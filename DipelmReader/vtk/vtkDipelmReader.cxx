// -------------------------------------------------------------------------- //
// UKRmol+ ParaView plugins                                                   //
// Jakub Benda (c) 2023                                                       //
// -------------------------------------------------------------------------- //

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>

// -------------------------------------------------------------------------- //

#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

// -------------------------------------------------------------------------- //

#include "vtkDipelmReader.h"

// -------------------------------------------------------------------------- //

vtkStandardNewMacro(vtkDipelmReader)

// -------------------------------------------------------------------------- //

const double Eh_eV = 27.21;

// -------------------------------------------------------------------------- //

vtkDipelmReader::vtkDipelmReader()
{
    SetNumberOfInputPorts(0);
}

vtkDipelmReader::~vtkDipelmReader()
{
}

void vtkDipelmReader::PrintSelf(std::ostream& os, vtkIndent indent)
{
    Superclass::PrintSelf(os, indent);

    os << indent << "File Name: " << (FileName ? FileName : "") << std::endl;
}

int vtkDipelmReader::CanReadFile(char const* name)
{
    return true;
}

int vtkDipelmReader::RequestData
(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector
)
{
    vtkInformation* outputInfo = outputVector->GetInformationObject(0);

    vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast
    (
        outputInfo->Get(vtkDataObject::DATA_OBJECT())
    );

    if (std::strlen(FileName) == 0)
        return 0;

    std::ifstream in(FileName, std::ios::binary);

    std::int64_t ionStateId = 0;
    std::int64_t neutralStateId = 0;
    std::int64_t nEnergies = 0;
    std::int64_t nThetas = 0;
    std::int64_t nPhis = 0;
    std::int64_t nComponents = 0;
    std::int64_t nAngles = 0;

    std::vector<double> energies, thetas, phis, dcs;

    std::streamsize intSize = sizeof(std::int64_t);
    std::streamsize floatSize = sizeof(double);

    in.read(reinterpret_cast<char*>(&ionStateId), intSize);
    in.read(reinterpret_cast<char*>(&neutralStateId), intSize);
    in.read(reinterpret_cast<char*>(&nEnergies), intSize);
    energies.resize(nEnergies);
    in.read(reinterpret_cast<char*>(&energies[0]), nEnergies*floatSize);
    in.read(reinterpret_cast<char*>(&nThetas), intSize);
    thetas.resize(nThetas);
    in.read(reinterpret_cast<char*>(&thetas[0]), nThetas*floatSize);
    in.read(reinterpret_cast<char*>(&nPhis), intSize);
    phis.resize(nPhis);
    in.read(reinterpret_cast<char*>(&phis[0]), nPhis*floatSize);
    in.read(reinterpret_cast<char*>(&nAngles), intSize);
    in.read(reinterpret_cast<char*>(&nComponents), intSize);
    in.read(reinterpret_cast<char*>(&nEnergies), intSize);
    dcs.resize(nAngles * nComponents * nEnergies);
    in.read(reinterpret_cast<char*>(&dcs[0]), dcs.size()*floatSize);

    vtkIdType nPoints = 0;
    vtkIdType nCells = 0;

    nPoints = nEnergies * (nPhis + 1) * nThetas;
    nCells = (nEnergies - 1) * nPhis * (nThetas - 1);

    std::vector<double> sinTheta(nThetas), cosTheta(nThetas), sinPhi(nPhis), cosPhi(nPhis);

    std::transform(thetas.begin(), thetas.end(), sinTheta.begin(), sin);
    std::transform(thetas.begin(), thetas.end(), cosTheta.begin(), cos);
    std::transform(phis.begin(), phis.end(), sinPhi.begin(), sin);
    std::transform(phis.begin(), phis.end(), cosPhi.begin(), cos);

    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(nPoints);

    vtkIdType iPoint = 0;

    for (vtkIdType iEnergy = 0; iEnergy < nEnergies; iEnergy++)
    {
        double energy = energies[iEnergy];

        if (EnergyUnit == 2)
            energy *= Eh_eV;

        energy += EnergyOffset;

        for (vtkIdType iPhi = 0; iPhi < nPhis; iPhi++)
        {
            for (vtkIdType iTheta = 0; iTheta < nThetas; iTheta++)
            {
                points->SetPoint
                (
                    iPoint++,
                    energy * sinTheta[iTheta] * cosPhi[iPhi],
                    energy * sinTheta[iTheta] * sinPhi[iPhi],
                    energy * cosTheta[iTheta]
                );
            }
        }
    }

    vtkNew<vtkCellArray> cellArray;
    vtkNew<vtkIdTypeArray> offsets, connectivity;

    offsets->SetNumberOfValues(nCells + 1);
    offsets->SetValue(0, 0);
    offsets->SetValue(nCells, 8*nCells);

    connectivity->SetNumberOfValues(8*nCells);

    vtkIdType iCell = 0;

    for (vtkIdType iEnergy = 0; iEnergy < nEnergies - 1; iEnergy++)
    {
        for (vtkIdType iPhi = 0; iPhi < nPhis; iPhi++)
        {
            for (vtkIdType iTheta = 0; iTheta < nThetas - 1; iTheta++)
            {
                offsets->SetValue(iCell, 8*iCell);

                connectivity->SetValue(8*iCell + 0, ((iEnergy + 0)*nPhis + (iPhi + 0) % nPhis)*nThetas + iTheta + 0);
                connectivity->SetValue(8*iCell + 1, ((iEnergy + 0)*nPhis + (iPhi + 0) % nPhis)*nThetas + iTheta + 1);
                connectivity->SetValue(8*iCell + 2, ((iEnergy + 0)*nPhis + (iPhi + 1) % nPhis)*nThetas + iTheta + 1);
                connectivity->SetValue(8*iCell + 3, ((iEnergy + 0)*nPhis + (iPhi + 1) % nPhis)*nThetas + iTheta + 0);
                connectivity->SetValue(8*iCell + 4, ((iEnergy + 1)*nPhis + (iPhi + 0) % nPhis)*nThetas + iTheta + 0);
                connectivity->SetValue(8*iCell + 5, ((iEnergy + 1)*nPhis + (iPhi + 0) % nPhis)*nThetas + iTheta + 1);
                connectivity->SetValue(8*iCell + 6, ((iEnergy + 1)*nPhis + (iPhi + 1) % nPhis)*nThetas + iTheta + 1);
                connectivity->SetValue(8*iCell + 7, ((iEnergy + 1)*nPhis + (iPhi + 1) % nPhis)*nThetas + iTheta + 0);

                iCell++;
            }
        }
    }

    cellArray->SetNumberOfCells(nCells);
    cellArray->SetData(offsets, connectivity);

    output->SetPoints(points);
    output->SetCells(VTK_HEXAHEDRON, cellArray);
    output->BuildLinks();

    for (int iComponent = 0; iComponent < nComponents; iComponent++)
    {
        const char* componentName[3] = { "Y pol", "Z pol", "X pol" };

        vtkFloatArray* array = vtkFloatArray::SafeDownCast(output->GetPointData()->GetAbstractArray(iComponent));

        if (array == nullptr)
        {
            vtkNew<vtkFloatArray> newArray;
            output->GetPointData()->AddArray(newArray);
            array = newArray;
        }

        array->SetNumberOfValues(nPoints);
        array->SetName(componentName[iComponent]);

        vtkIdType iPoint = 0;

        for (vtkIdType iEnergy = 0; iEnergy < nEnergies; iEnergy++)
        {
            for (vtkIdType iPhi = 0; iPhi < nPhis; iPhi++)
            {
                for (vtkIdType iTheta = 0; iTheta < nThetas; iTheta++)
                {
                    array->SetValue(iPoint++, dcs[((iEnergy*nComponents + iComponent)*nPhis + iPhi)*nThetas + iTheta]);
                }
            }
        }
    }

    return 1;
}

// -------------------------------------------------------------------------- //
