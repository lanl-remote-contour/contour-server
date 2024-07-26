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
#include "vtkContourPreFilter.h"

#include <vtkDataArraySelection.h>
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkZLibDataCompressor.h>

#include <filesystem>
#include <string>

size_t CompressMap(const std::unordered_map<int, float>& map)
{
  vtkNew<vtkZLibDataCompressor> compressor;
  size_t sz0 = map.size() << 3;
  auto* const buf0 = new unsigned char[sz0];
  auto* const buf1 = new unsigned char[sz0];
  unsigned char* ptr = buf0;
  for (auto const& [key, val] : map)
  {
    memcpy(ptr, (unsigned char*)&key, 4);
    ptr += 4;
    memcpy(ptr, (unsigned char*)&val, 4);
    ptr += 4;
  }
  size_t sz1 = compressor->Compress(buf0, sz0, buf1, sz0);
  delete[] buf0;
  delete[] buf1;
  return sz1;
}

void ProcessContourValue(vtkContourPreFilter* pc, double value)
{
  pc->SetValue(0, value);
  pc->Update();
  std::cout << value << " " << pc->GetResult().size() << " " << CompressMap(pc->GetResult())
            << std::endl;
}

void ProcessFile(const char* fileName, const char* arrayName)
{
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(fileName);
  reader->UpdateInformation();
  vtkDataArraySelection* das = reader->GetPointDataArraySelection();
  das->DisableAllArrays();
  das->EnableArray(arrayName);

  vtkNew<vtkContourPreFilter> pc;
  pc->SetInputConnection(reader->GetOutputPort());
  pc->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FieldAssociations::FIELD_ASSOCIATION_POINTS, arrayName);

  for (int i = 1; i <= 9; i++)
    ProcessContourValue(pc, i / 10.0);
}

int main(int argc, char* argv[])
{
  char arrayName[] = "v02";
  for (const auto& entry : std::filesystem::directory_iterator(argv[1]))
  {
    std::cout << entry.path() << std::endl;
    ProcessFile(entry.path().c_str(), arrayName);
  }

  return 0;
}
