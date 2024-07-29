/*
 * Copyright (c) 2024 Triad National Security, LLC, as operator of Los Alamos
 * National Laboratory with the U.S. Department of Energy/National Nuclear
 * Security Administration. All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * with the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of TRIAD, Los Alamos National Laboratory, LANL, the
 *    U.S. Government, nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Copyright (c) 2024 Triad National Security, LLC, as operator of Los Alamos
 * National Laboratory with the U.S. Department of Energy/National Nuclear
 * Security Administration. All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * with the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of TRIAD, Los Alamos National Laboratory, LANL, the
 *    U.S. Government, nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "vtkPostFilteringSynchronizedTemplates3D.h"

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>

#include <rpc/client.h>

#include <chrono>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>

void MakeContour(
  rpc::client* cli, const std::string& fileName, const std::string& arrayName, double value)
{
  vtkNew<vtkFloatArray> scalars;
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(500 * 500 * 500);
  scalars->FillValue(-1);

  vtkNew<vtkImageData> image;
  image->SetOrigin(-2300000, -500000, -1200000);
  image->SetSpacing(9218.4368737, 5611.2224449, 4809.6192385);
  image->SetExtent(0, 499, 0, 499, 0, 499);
  image->GetPointData()->SetScalars(scalars);

  auto t0 = std::chrono::high_resolution_clock::now();

  std::unordered_map<int, float> map =
    cli->call("LoadContour", fileName, arrayName, value).as<std::unordered_map<int, float>>();

  auto t1 = std::chrono::high_resolution_clock::now();

  float* ptr = scalars->GetPointer(0);
  for (auto const& [key, val] : map)
  {
    ptr[key] = val;
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  vtkNew<vtkPostFilteringSynchronizedTemplates3D> contour;
  contour->SetInputData(image);
  contour->SetValue(0, value);
  contour->Update();

  auto t3 = std::chrono::high_resolution_clock::now();

  std::cout << "time: " << std::chrono::duration<double>(t1 - t0).count() << " "
            << std::chrono::duration<double>(t2 - t1).count() << " "
            << std::chrono::duration<double>(t3 - t2).count() << std::endl;
}

int main(int argc, char* argv[])
{
  const char* arr = "v03";
  double value = 0.9;
  int c;
  while ((c = getopt(argc, argv, "a:c:")) != -1)
  {
    switch (c)
    {
      case 'a':
        arr = optarg;
        break;
      case 'c':
        value = atof(optarg);
        break;
      default:
        std::cerr << "Use -a to specify array name and -c to specify contour value" << std::endl;
        exit(EXIT_FAILURE);
    }
  }
  argc -= optind;
  argv += optind;
  if (argc < 2)
  {
    std::cerr << "Lack server_ip and file_name" << std::endl;
    exit(EXIT_FAILURE);
  }
  rpc::client cli(argv[0], 8080);
  std::cout << "vtk file: " << argv[1] << std::endl;
  std::cout << "contour value: " << value << std::endl;
  std::cout << "array: " << arr << std::endl;
  std::cout << "Connecting to " << argv[0] << "..." << std::endl;
  while (cli.get_connection_state() != rpc::client::connection_state::connected)
  {
    // Wait
  }
  std::cout << "Done!" << std::endl;
  MakeContour(&cli, argv[1], arr, value);
  return 0;
}
