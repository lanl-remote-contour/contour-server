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

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdListCollection.h>
#include <vtkInformation.h>
#include <vtkInformationIntegerVectorKey.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkShortArray.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStructuredPoints.h>
#include <vtkSynchronizedTemplates3D.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkUnsignedShortArray.h>

#include <cmath>

VTK_ABI_NAMESPACE_BEGIN
vtkStandardNewMacro(vtkPostFilteringSynchronizedTemplates3D);

//------------------------------------------------------------------------------
// Description:
// Construct object with initial scalar range (0,1) and single contour value of 0.0.
vtkPostFilteringSynchronizedTemplates3D::vtkPostFilteringSynchronizedTemplates3D()
{
  this->ContourValues = vtkContourValues::New();
  this->ArrayComponent = 0;

  // by default process active point scalars
  this->SetInputArrayToProcess(
    0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
}

//------------------------------------------------------------------------------
vtkPostFilteringSynchronizedTemplates3D::~vtkPostFilteringSynchronizedTemplates3D()
{
  this->ContourValues->Delete();
}

//------------------------------------------------------------------------------
// Overload standard modified time function. If contour values are modified,
// then this object is modified as well.
vtkMTimeType vtkPostFilteringSynchronizedTemplates3D::GetMTime()
{
  vtkMTimeType mTime = this->Superclass::GetMTime();
  vtkMTimeType mTime2 = this->ContourValues->GetMTime();

  mTime = (mTime2 > mTime ? mTime2 : mTime);
  return mTime;
}

//------------------------------------------------------------------------------
static void vtkPostFilteringSynchronizedTemplates3DInitializeOutput(
  int* ext, vtkImageData* input, vtkPolyData* o, vtkFloatArray* scalars, vtkDataArray* inScalars)
{
  vtkPoints* newPts;
  vtkCellArray* newPolys;

  const vtkIdType numCells = static_cast<vtkIdType>(ext[1] - ext[0] + 1) *
    static_cast<vtkIdType>(ext[3] - ext[2] + 1) * static_cast<vtkIdType>(ext[5] - ext[4] + 1);

  vtkIdType estimatedSize = (vtkIdType)pow(static_cast<double>(numCells), .75);
  if (estimatedSize < 1024)
  {
    estimatedSize = 1024;
  }
  newPts = vtkPoints::New();
  newPts->Allocate(estimatedSize, estimatedSize);
  newPolys = vtkCellArray::New();
  newPolys->AllocateEstimate(estimatedSize, 3);

  o->GetPointData()->CopyAllOn();
  // It is more efficient to just create the scalar array
  // rather than redundantly interpolate the scalars.
  if (input->GetPointData()->GetScalars() == inScalars)
  {
    o->GetPointData()->CopyScalarsOff();
  }
  else
  {
    o->GetPointData()->CopyFieldOff(inScalars->GetName());
  }

  if (scalars)
  {
    // A temporary name.
    scalars->SetName("Scalars");
  }

  o->GetPointData()->InterpolateAllocate(input->GetPointData(), estimatedSize, estimatedSize / 2);
  o->GetCellData()->CopyAllocate(input->GetCellData(), estimatedSize, estimatedSize / 2);

  o->SetPoints(newPts);
  newPts->Delete();

  o->SetPolys(newPolys);
  newPolys->Delete();
}

//------------------------------------------------------------------------------
//
// Contouring filter specialized for images
//
template <class T>
void ContourImage(vtkPostFilteringSynchronizedTemplates3D* self, int* exExt, vtkImageData* data,
  vtkPolyData* output, T* ptr, vtkDataArray* inScalars)
{
  int* inExt = data->GetExtent();
  vtkIdType xdim = exExt[1] - exExt[0] + 1;
  vtkIdType ydim = exExt[3] - exExt[2] + 1;
  double* values = self->GetValues();
  vtkIdType numContours = self->GetNumberOfContours();
  T *inPtrX, *inPtrY, *inPtrZ;
  T *s0, *s1, *s2, *s3;
  int xMin, xMax, yMin, yMax, zMin, zMax;
  vtkIdType xInc, yInc, zInc;
  double* origin = data->GetOrigin();
  double* spacing = data->GetSpacing();
  vtkIdType *isect1Ptr, *isect2Ptr;
  double y, z, t;
  int i, j, k;
  vtkIdType zstep, yisectstep;
  vtkIdType offsets[12];
  vtkTypeBool ComputeScalars = self->GetComputeScalars();
  int* tablePtr;
  vtkIdType idx;
  int vidx;
  double x[3], xz[3];
  vtkIdType v0, v1, v2, v3;
  vtkIdType ptIds[3];
  double value;
  // We need to know the edgePointId's for interpolating attributes.
  vtkIdType edgePtId, inCellId, outCellId;
  vtkPointData* inPD = data->GetPointData();
  vtkCellData* inCD = data->GetCellData();
  vtkPointData* outPD = output->GetPointData();
  vtkCellData* outCD = output->GetCellData();
  vtkFloatArray* newScalars = nullptr;
  vtkPoints* newPts;
  vtkCellArray* newPolys;
  ptr += self->GetArrayComponent();
  bool abort = false;

  if (ComputeScalars)
  {
    newScalars = vtkFloatArray::New();
  }
  vtkPostFilteringSynchronizedTemplates3DInitializeOutput(
    exExt, data, output, newScalars, inScalars);
  newPts = output->GetPoints();
  newPolys = output->GetPolys();

  // this is an exploded execute extent.
  xMin = exExt[0];
  xMax = exExt[1];
  yMin = exExt[2];
  yMax = exExt[3];
  zMin = exExt[4];
  zMax = exExt[5];

  // increments to move through scalars. Compute these ourself because
  // we may be contouring an array other than scalars.
  xInc = inScalars->GetNumberOfComponents();
  yInc = xInc * (inExt[1] - inExt[0] + 1);
  zInc = yInc * (inExt[3] - inExt[2] + 1);

  // Kens increments, probably to do with edge array
  zstep = xdim * ydim;
  yisectstep = xdim * 3;
  // compute offsets probably how to get to the edges in the edge array.
  offsets[0] = -xdim * 3;
  offsets[1] = -xdim * 3 + 1;
  offsets[2] = -xdim * 3 + 2;
  offsets[3] = -xdim * 3 + 4;
  offsets[4] = -xdim * 3 + 5;
  offsets[5] = 0;
  offsets[6] = 2;
  offsets[7] = 5;
  offsets[8] = (zstep - xdim) * 3;
  offsets[9] = (zstep - xdim) * 3 + 1;
  offsets[10] = (zstep - xdim) * 3 + 4;
  offsets[11] = zstep * 3;

  // allocate storage array
  vtkIdType* isect1 = new vtkIdType[xdim * ydim * 3 * 2];
  // set impossible edges to -1
  for (i = 0; i < ydim; i++)
  {
    isect1[(i + 1) * xdim * 3 - 3] = -1;
    isect1[(i + 1) * xdim * 3 * 2 - 3] = -1;
  }
  for (i = 0; i < xdim; i++)
  {
    isect1[((ydim - 1) * xdim + i) * 3 + 1] = -1;
    isect1[((ydim - 1) * xdim + i) * 3 * 2 + 1] = -1;
  }

  int checkAbortInterval = std::min((zMax - zMin) / 10 + 1, 1000);

  // for each contour
  for (vidx = 0; vidx < numContours && !abort; vidx++)
  {
    value = values[vidx];
    inPtrZ = ptr;

    //==================================================================
    for (k = zMin; k <= zMax; k++)
    {
      self->UpdateProgress(
        (double)vidx / numContours + (k - zMin) / ((zMax - zMin + 1.0) * numContours));
      if (k % checkAbortInterval == 0 && self->CheckAbort())
      {
        abort = true;
        break;
      }
      z = origin[2] + spacing[2] * k;
      x[2] = z;

      // swap the buffers
      if (k % 2)
      {
        offsets[8] = (zstep - xdim) * 3;
        offsets[9] = (zstep - xdim) * 3 + 1;
        offsets[10] = (zstep - xdim) * 3 + 4;
        offsets[11] = zstep * 3;
        isect1Ptr = isect1;
        isect2Ptr = isect1 + xdim * ydim * 3;
      }
      else
      {
        offsets[8] = (-zstep - xdim) * 3;
        offsets[9] = (-zstep - xdim) * 3 + 1;
        offsets[10] = (-zstep - xdim) * 3 + 4;
        offsets[11] = -zstep * 3;
        isect1Ptr = isect1 + xdim * ydim * 3;
        isect2Ptr = isect1;
      }

      inPtrY = inPtrZ;
      for (j = yMin; j <= yMax; j++)
      {
        // Should not impact performance here/
        edgePtId = (xMin - inExt[0]) * xInc + (j - inExt[2]) * yInc + (k - inExt[4]) * zInc;
        // Increments are different for cells.  Since the cells are not
        // contoured until the second row of templates, subtract 1 from
        // i,j,and k.  Note: first cube is formed when i=0, j=1, and k=1.
        inCellId = (xMin - inExt[0]) +
          (inExt[1] - inExt[0]) * ((j - inExt[2] - 1) + (k - inExt[4] - 1) * (inExt[3] - inExt[2]));

        y = origin[1] + j * spacing[1];
        xz[1] = y;
        s1 = inPtrY;
        v1 = (*s1 < value ? 0 : 1);

        inPtrX = inPtrY;
        for (i = xMin; i <= xMax; i++)
        {
          s0 = s1;
          v0 = v1;
          *isect2Ptr = -1;
          *(isect2Ptr + 1) = -1;
          *(isect2Ptr + 2) = -1;
          if (i < xMax)
          {
            s1 = (inPtrX + xInc);
            v1 = (*s1 < value ? 0 : 1);
            if ((*s0 != -1) && (*s1 != -1) && (v0 ^ v1))
            {
              // watch for degenerate points
              if (*s0 == value)
              {
                if (i > xMin && *(isect2Ptr - 3) > -1)
                {
                  *isect2Ptr = *(isect2Ptr - 3);
                }
                else if (j > yMin && *(isect2Ptr - yisectstep + 1) > -1)
                {
                  *isect2Ptr = *(isect2Ptr - yisectstep + 1);
                }
                else if (k > zMin && *(isect1Ptr + 2) > -1)
                {
                  *isect2Ptr = *(isect1Ptr + 2);
                }
              }
              else if (*s1 == value)
              {
                if (j > yMin && *(isect2Ptr - yisectstep + 4) > -1)
                {
                  *isect2Ptr = *(isect2Ptr - yisectstep + 4);
                }
                else if (k > zMin && i < xMax && *(isect1Ptr + 5) > -1)
                {
                  *isect2Ptr = *(isect1Ptr + 5);
                }
              }
              // if the edge has not been set yet then it is a new point
              if (*isect2Ptr == -1)
              {
                t = (value - (double)(*s0)) / ((double)(*s1) - (double)(*s0));
                x[0] = origin[0] + spacing[0] * (i + t);
                x[1] = y;
                *isect2Ptr = newPts->InsertNextPoint(x);
                if (newScalars)
                  newScalars->InsertNextTuple(&value);
                outPD->InterpolateEdge(inPD, *isect2Ptr, edgePtId, edgePtId + 1, t);
              }
            }
          }
          if (j < yMax)
          {
            s2 = (inPtrX + yInc);
            v2 = (*s2 < value ? 0 : 1);
            if ((*s0 != -1) && (*s2 != -1) && (v0 ^ v2))
            {
              if (*s0 == value)
              {
                if (*isect2Ptr > -1)
                {
                  *(isect2Ptr + 1) = *isect2Ptr;
                }
                else if (i > xMin && *(isect2Ptr - 3) > -1)
                {
                  *(isect2Ptr + 1) = *(isect2Ptr - 3);
                }
                else if (j > yMin && *(isect2Ptr - yisectstep + 1) > -1)
                {
                  *(isect2Ptr + 1) = *(isect2Ptr - yisectstep + 1);
                }
                else if (k > zMin && *(isect1Ptr + 2) > -1)
                {
                  *(isect2Ptr + 1) = *(isect1Ptr + 2);
                }
              }
              else if (*s2 == value && k > zMin && *(isect1Ptr + yisectstep + 2) > -1)
              {
                *(isect2Ptr + 1) = *(isect1Ptr + yisectstep + 2);
              }
              // if the edge has not been set yet then it is a new point
              if (*(isect2Ptr + 1) == -1)
              {
                t = (value - (double)(*s0)) / ((double)(*s2) - (double)(*s0));
                x[0] = origin[0] + spacing[0] * i;
                x[1] = y + spacing[1] * t;
                *(isect2Ptr + 1) = newPts->InsertNextPoint(x);
                if (newScalars)
                  newScalars->InsertNextTuple(&value);
                outPD->InterpolateEdge(inPD, *(isect2Ptr + 1), edgePtId, edgePtId + yInc, t);
              }
            }
          }
          if (k < zMax)
          {
            s3 = (inPtrX + zInc);
            v3 = (*s3 < value ? 0 : 1);
            if ((*s0 != -1) && (*s3 != -1) && (v0 ^ v3))
            {
              if (*s0 == value)
              {
                if (*isect2Ptr > -1)
                {
                  *(isect2Ptr + 2) = *isect2Ptr;
                }
                else if (*(isect2Ptr + 1) > -1)
                {
                  *(isect2Ptr + 2) = *(isect2Ptr + 1);
                }
                else if (i > xMin && *(isect2Ptr - 3) > -1)
                {
                  *(isect2Ptr + 2) = *(isect2Ptr - 3);
                }
                else if (j > yMin && *(isect2Ptr - yisectstep + 1) > -1)
                {
                  *(isect2Ptr + 2) = *(isect2Ptr - yisectstep + 1);
                }
                else if (k > zMin && *(isect1Ptr + 2) > -1)
                {
                  *(isect2Ptr + 2) = *(isect1Ptr + 2);
                }
              }
              if (*(isect2Ptr + 2) == -1)
              {
                t = (value - (double)(*s0)) / ((double)(*s3) - (double)(*s0));
                xz[0] = origin[0] + spacing[0] * i;
                xz[2] = z + spacing[2] * t;
                *(isect2Ptr + 2) = newPts->InsertNextPoint(xz);
                if (newScalars)
                  newScalars->InsertNextTuple(&value);
                outPD->InterpolateEdge(inPD, *(isect2Ptr + 2), edgePtId, edgePtId + zInc, t);
              }
            }
          }
          // To keep track of ids for interpolating attributes.
          ++edgePtId;

          // now add any polys that need to be added
          // basically look at the isect values,
          // form an index and lookup the polys
          if ((*s0 != -1) && j > yMin && i < xMax && k > zMin)
          {
            idx = (v0 ? 4096 : 0);
            idx = idx + (*(isect1Ptr - yisectstep) > -1 ? 2048 : 0);
            idx = idx + (*(isect1Ptr - yisectstep + 1) > -1 ? 1024 : 0);
            idx = idx + (*(isect1Ptr - yisectstep + 2) > -1 ? 512 : 0);
            idx = idx + (*(isect1Ptr - yisectstep + 4) > -1 ? 256 : 0);
            idx = idx + (*(isect1Ptr - yisectstep + 5) > -1 ? 128 : 0);
            idx = idx + (*(isect1Ptr) > -1 ? 64 : 0);
            idx = idx + (*(isect1Ptr + 2) > -1 ? 32 : 0);
            idx = idx + (*(isect1Ptr + 5) > -1 ? 16 : 0);
            idx = idx + (*(isect2Ptr - yisectstep) > -1 ? 8 : 0);
            idx = idx + (*(isect2Ptr - yisectstep + 1) > -1 ? 4 : 0);
            idx = idx + (*(isect2Ptr - yisectstep + 4) > -1 ? 2 : 0);
            idx = idx + (*(isect2Ptr) > -1 ? 1 : 0);

            tablePtr =
              VTK_SYNCHRONIZED_TEMPLATES_3D_TABLE_2 + VTK_SYNCHRONIZED_TEMPLATES_3D_TABLE_1[idx];

            while (*tablePtr != -1)
            {
              ptIds[0] = *(isect1Ptr + offsets[*tablePtr]);
              tablePtr++;
              ptIds[1] = *(isect1Ptr + offsets[*tablePtr]);
              tablePtr++;
              ptIds[2] = *(isect1Ptr + offsets[*tablePtr]);
              tablePtr++;
              if (ptIds[0] != ptIds[1] && ptIds[0] != ptIds[2] && ptIds[1] != ptIds[2])
              {
                outCellId = newPolys->InsertNextCell(3, ptIds);
                outCD->CopyData(inCD, inCellId, outCellId);
              }
            }
          }
          inPtrX += xInc;
          isect2Ptr += 3;
          isect1Ptr += 3;
          // To keep track of ids for copying cell attributes..
          ++inCellId;
        }
        inPtrY += yInc;
      }
      inPtrZ += zInc;
    }
  }
  delete[] isect1;

  if (newScalars)
  {
    // Lets set the name of the scalars here.
    if (inScalars)
    {
      newScalars->SetName(inScalars->GetName());
    }
    idx = output->GetPointData()->AddArray(newScalars);
    output->GetPointData()->SetActiveAttribute(idx, vtkDataSetAttributes::SCALARS);
    newScalars->Delete();
    newScalars = nullptr;
  }
}

//------------------------------------------------------------------------------
//
// Contouring filter specialized for images (or slices from images)
//
void vtkPostFilteringSynchronizedTemplates3D::ThreadedExecute(
  vtkImageData* data, vtkInformation* inInfo, vtkInformation* outInfo, vtkDataArray* inScalars)
{
  void* ptr;
  vtkPolyData* output;

  vtkDebugMacro(<< "Executing 3D structured contour");

  output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

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
  switch (inScalars->GetDataType())
  {
    vtkTemplateMacro(ContourImage(this, exExt, data, output, (VTK_TT*)ptr, inScalars));
  }
}

//------------------------------------------------------------------------------
int vtkPostFilteringSynchronizedTemplates3D::RequestData(
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkImageData* input = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // to be safe recompute the
  this->RequestUpdateExtent(request, inputVector, outputVector);

  vtkDataArray* inScalars = this->GetInputArrayToProcess(0, inputVector);

  // Just call the threaded execute directly.
  this->ThreadedExecute(input, inInfo, outInfo, inScalars);

  output->Squeeze();

  return 1;
}

//------------------------------------------------------------------------------
int vtkPostFilteringSynchronizedTemplates3D::FillInputPortInformation(int, vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//------------------------------------------------------------------------------
void vtkPostFilteringSynchronizedTemplates3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  this->ContourValues->PrintSelf(os, indent.GetNextIndent());
  os << indent << "Compute Scalars: " << (this->ComputeScalars ? "On\n" : "Off\n");
  os << indent << "ArrayComponent: " << this->ArrayComponent << endl;
}

VTK_ABI_NAMESPACE_END
