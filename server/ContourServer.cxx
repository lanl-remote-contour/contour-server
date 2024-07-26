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

#include <rpc/server.h>

std::unordered_map<int, float> LoadContour(
  const std::string& fileName, const std::string& arrayName, double value)
{
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(fileName.c_str());
  reader->UpdateInformation();
  vtkDataArraySelection* das = reader->GetPointDataArraySelection();
  das->DisableAllArrays();
  das->EnableArray(arrayName.c_str());

  vtkNew<vtkContourPreFilter> pc;
  pc->SetInputConnection(reader->GetOutputPort());
  pc->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FieldAssociations::FIELD_ASSOCIATION_POINTS, arrayName.c_str());
  pc->SetValue(0, value);
  pc->Update();

  std::unordered_map<int, float> result;
  result.swap(pc->GetResult());
  return result;
}

int main(int argc, char* argv[])
{
  rpc::server srv(8080);
  srv.bind("LoadContour", LoadContour);
  srv.run();
  return 0;
}
