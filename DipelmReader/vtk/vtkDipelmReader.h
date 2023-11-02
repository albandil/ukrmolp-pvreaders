// -------------------------------------------------------------------------- //
// UKRmol+ ParaView plugins                                                   //
// Jakub Benda (c) 2023                                                       //
// -------------------------------------------------------------------------- //

#ifndef __vtkDipelmReader_h__
#define __vtkDipelmReader_h__

// -------------------------------------------------------------------------- //

#include <vtkUnstructuredGridAlgorithm.h>

// -------------------------------------------------------------------------- //

class VTK_EXPORT vtkDipelmReader : public vtkUnstructuredGridAlgorithm
{
    public:

        static vtkDipelmReader* New();
        vtkTypeMacro(vtkDipelmReader, vtkAlgorithm);
        void PrintSelf(std::ostream& os, vtkIndent indent) override;

        virtual int CanReadFile(char const* name);

        vtkSetMacro(EnergyUnit, int);
        vtkGetMacro(EnergyUnit, int);

        vtkSetMacro(EnergyOffset, double);
        vtkGetMacro(EnergyOffset, double);

        vtkGetStringMacro(FileName)
        vtkSetStringMacro(FileName)

    protected:

        vtkDipelmReader();
        ~vtkDipelmReader();

        int RequestData
        (
            vtkInformation *request,
            vtkInformationVector **inputVector,
            vtkInformationVector *outputVector
        ) override;

    private:

        char *FileName{};
        int EnergyUnit{};
        double EnergyOffset{};
};

// -------------------------------------------------------------------------- //

#endif

// -------------------------------------------------------------------------- //
