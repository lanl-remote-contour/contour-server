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

#include <vtkContourFilter.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>

#include <chrono>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>

void LoadContour(const char* fileName, const char* arrayName, double value)
{
  auto t0 = std::chrono::high_resolution_clock::now();

  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(fileName);
  reader->UpdateInformation();
  vtkDataArraySelection* das = reader->GetPointDataArraySelection();
  das->DisableAllArrays();
  das->EnableArray(arrayName);
  reader->Update();

  auto t1 = std::chrono::high_resolution_clock::now();

  vtkNew<vtkContourFilter> contour;
  contour->SetInputConnection(reader->GetOutputPort());
  contour->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FieldAssociations::FIELD_ASSOCIATION_POINTS, arrayName);
  contour->SetValue(0, value);
  contour->Update();

  auto t2 = std::chrono::high_resolution_clock::now();

  std::cout << "time: " << std::chrono::duration<double>(t1 - t0).count() << " "
            << std::chrono::duration<double>(t2 - t1).count() << std::endl;
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
  if (!argc)
  {
    std::cerr << "Lack target vtk filename" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::cout << "vtk file: " << argv[1] << std::endl;
  std::cout << "contour value: " << value << std::endl;
  std::cout << "array: " << arr << std::endl;
  LoadContour(argv[1], arr, value);
  return 0;
}
