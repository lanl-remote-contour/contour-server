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

#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDataSetAttributes.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkContourPreFilter);

//------------------------------------------------------------------------------
vtkContourPreFilter::vtkContourPreFilter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(0);
  this->ContourValues = vtkContourValues::New();
  this->ArrayComponent = 0;

  this->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
}

//------------------------------------------------------------------------------
vtkContourPreFilter::~vtkContourPreFilter()
{
  this->ContourValues->Delete();
}

//------------------------------------------------------------------------------
vtkMTimeType vtkContourPreFilter::GetMTime()
{
  vtkMTimeType mTime = this->Superclass::GetMTime();
  vtkMTimeType mTime2 = this->ContourValues->GetMTime();

  mTime = (mTime2 > mTime ? mTime2 : mTime);
  return mTime;
}

template <class T>
void PreProcessImage(
  vtkContourPreFilter* self, int* exExt, vtkImageData* data, T* ptr, vtkDataArray* inScalars)
{
  std::unordered_map<int, T> result;
  int* inExt = data->GetExtent();
  ptr += self->GetArrayComponent();
  T *inPtrX, *inPtrY, *inPtrZ;
  vtkIdType v0, v1, v2, v3;
  T *s0, *s1, *s2, *s3;
  int xMin, xMax, yMin, yMax, zMin, zMax;
  vtkIdType xInc, yInc, zInc;
  int i, j, k;

  xMin = exExt[0];
  xMax = exExt[1];
  yMin = exExt[2];
  yMax = exExt[3];
  zMin = exExt[4];
  zMax = exExt[5];

  xInc = inScalars->GetNumberOfComponents();
  yInc = xInc * (inExt[1] - inExt[0] + 1);
  zInc = yInc * (inExt[3] - inExt[2] + 1);

  double* values = self->GetValues();
  const double value = values[0];

  inPtrZ = ptr;
  for (k = zMin; k <= zMax; k++)
  {
    inPtrY = inPtrZ;
    for (j = yMin; j <= yMax; j++)
    {
      s1 = inPtrY;
      v1 = (*s1 < value ? 0 : 1);

      inPtrX = inPtrY;
      for (i = xMin; i <= xMax; i++)
      {
        s0 = s1;
        v0 = v1;

        if (i < xMax)
        {
          s1 = (inPtrX + xInc);
          v1 = (*s1 < value ? 0 : 1);
          if (v0 ^ v1)
          {
            result[(inPtrX - ptr) / xInc] = *s0;
            result[(inPtrX + xInc - ptr) / xInc] = *s1;
          }
        }

        if (j < yMax)
        {
          s2 = (inPtrX + yInc);
          v2 = (*s2 < value ? 0 : 1);

          if (v0 ^ v2)
          {
            result[(inPtrX - ptr) / xInc] = *s0;
            result[(inPtrX + yInc - ptr) / xInc] = *s2;
          }
        }

        if (k < zMax)
        {
          s3 = (inPtrX + zInc);
          v3 = (*s3 < value ? 0 : 1);

          if (v0 ^ v3)
          {
            result[(inPtrX - ptr) / xInc] = *s0;
            result[(inPtrX + zInc - ptr) / xInc] = *s3;
          }
        }

        inPtrX += xInc;
      }
      inPtrY += yInc;
    }
    inPtrZ += zInc;
  }

  result.swap(self->GetResult());
}

//------------------------------------------------------------------------------
vtkTypeBool vtkContourPreFilter::ProcessRequest(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
  {
    return this->RequestData(request, inputVector, outputVector);
  }

  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_TIME()))
  {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
  {
    return 1;
  }

  if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
  {
    return 1;
  }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//------------------------------------------------------------------------------
void vtkContourPreFilter::ThreadedExecute(
  vtkImageData* data, vtkInformation* inInfo, vtkDataArray* inScalars)
{
  void* ptr;

  vtkDebugMacro(<< "Executing contour prefilter");

  int* inExt = data->GetExtent();
  int exExt[6];
  inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), exExt);

  for (int i = 0; i < 3; i++)
  {
    if (inExt[2 * i] > exExt[2 * i])
    {
      exExt[2 * i] = inExt[2 * i];
    }
    if (inExt[2 * i + 1] < exExt[2 * i + 1])
    {
      exExt[2 * i + 1] = inExt[2 * i + 1];
    }
  }
  if (exExt[0] >= exExt[1] || exExt[2] >= exExt[3] || exExt[4] >= exExt[5])
  {
    vtkDebugMacro(<< "3D structured contours requires 3D data");
    return;
  }

  //
  // Check data type and execute appropriate function
  //
  if (inScalars == nullptr)
  {
    vtkDebugMacro("No scalars for contouring.");
    return;
  }
  int numComps = inScalars->GetNumberOfComponents();

  if (this->ArrayComponent >= numComps)
  {
    vtkErrorMacro("Scalars have " << numComps
                                  << " components. "
                                     "ArrayComponent must be smaller than "
                                  << numComps);

    return;
  }

  ptr = data->GetArrayPointerForExtent(inScalars, exExt);
  PreProcessImage(this, exExt, data, (float*)ptr, inScalars);
}

//------------------------------------------------------------------------------
int vtkContourPreFilter::RequestUpdateExtent(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector*)
{
  vtkInformation* inputInfo = inputVector[0]->GetInformationObject(0);
  inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}

//------------------------------------------------------------------------------
int vtkContourPreFilter::RequestData(
  vtkInformation*, vtkInformationVector** inputVector, vtkInformationVector*)
{
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray* inScalars = this->GetInputArrayToProcess(0, inputVector);
  this->ThreadedExecute(input, inInfo, inScalars);

  return 1;
}

//------------------------------------------------------------------------------
int vtkContourPreFilter::FillInputPortInformation(int, vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//------------------------------------------------------------------------------
void vtkContourPreFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  this->ContourValues->PrintSelf(os, indent.GetNextIndent());
  os << indent << "ArrayComponent: " << this->ArrayComponent << endl;
}

VTK_ABI_NAMESPACE_END
