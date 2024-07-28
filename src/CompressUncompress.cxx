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

#include <vtkNew.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>

#include <filesystem>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

void ConvertTo(const char* input, vtkXMLImageDataWriter::CompressorType type, const char* ext)
{
  vtkNew<vtkXMLImageDataReader> reader;
  reader->SetFileName(input);
  vtkNew<vtkXMLImageDataWriter> writer;
  writer->SetInputConnection(reader->GetOutputPort());
  writer->SetCompressorType(type);
  writer->SetFileName(std::filesystem::path(input).replace_extension(ext).c_str());
  writer->Write();
}

int main(int argc, char* argv[])
{
  vtkXMLImageDataWriter::CompressorType type = vtkXMLImageDataWriter::CompressorType::NONE;
  const char* ext = "raw";
  int c;
  while ((c = getopt(argc, argv, "glr")) != -1)
  {
    switch (c)
    {
      case 'g':
        type = vtkXMLImageDataWriter::CompressorType::ZLIB;
        ext = "gz";
        break;
      case 'l':
        type = vtkXMLImageDataWriter::CompressorType::LZ4;
        ext = "lz4";
        break;
      case 'r':
        type = vtkXMLImageDataWriter::CompressorType::NONE;
        ext = "raw";
        break;
      default:
        fprintf(stderr, "Use -g (ZLib), -l (LZ4), or -r (RAW) to specify compression types\n");
        exit(EXIT_FAILURE);
    }
  }
  argc -= optind;
  argv += optind;
  if (!argc)
  {
    fprintf(stderr, "Lack input filename\n");
    exit(EXIT_FAILURE);
  }
  ConvertTo(argv[0], type, ext);
  return 0;
}
