#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

void selectEvents(const char* inputFile, const char* outputFilename, int nmax) {
    // Open the input file
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        printf("Error: Unable to open input file %s\n", inputFile);
        return;
    }

    // Get the input tree
    TTree *inputTree = (TTree*)file->Get("Truth_Info");
    if (!inputTree) {
        printf("Error: Unable to retrieve input Truth_Info\n");
        file->Close();
        return;
    }

    // Create a new file
    TFile *outputFile = new TFile(outputFilename, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        printf("Error: Unable to create output file %s\n", outputFilename);
        file->Close();
        return;
    }

    // Clone the input tree
    TTree *outputTree = inputTree->CloneTree(0); // Create an empty tree with the same branches

    // Define the cut condition
    // Example: Select events where positionX, positionY, and positionZ are within a certain range
    const double LAr_Start_Exact[] = {-3478.48, -2166.71, 4179.24};
    const double LAr_End_Exact[] = {3478.48, 829.282, 9135.88};
    Double_t minX = LAr_Start_Exact[0], maxX = LAr_End_Exact[0];
    Double_t minY = LAr_Start_Exact[1], maxY = LAr_End_Exact[1];
    Double_t minZ = LAr_Start_Exact[2], maxZ = LAr_End_Exact[2];

    // Set up TTreeReader for input tree
    //TTreeReader reader(inputTree);
    //inputTree->SetBranchStatus("*", true);
    //TTreeReaderArray<float> X4(reader, "LeptonX4");
    
    float X4[4];
    inputTree->SetBranchAddress("LeptonX4", &X4);
    
    // Loop over events
    long n = 0;
    while (n < inputTree->GetEntries() && (nmax < 0 || n < nmax)) {
        inputTree->GetEntry(n);
        // Access the position components
        Double_t posX = X4[0];
        Double_t posY = X4[1];
        Double_t posZ = X4[2];
        // Apply the cut condition
        if (posX >= minX && posX <= maxX &&
            posY >= minY && posY <= maxY &&
            posZ >= minZ && posZ <= maxZ) {
            outputTree->Fill(); // Fill the output tree with selected events
        }
        n += 1;
    }

    // Write the output tree to the output file
    outputTree->Write();
    
    // Clean up
    outputFile->Close();
    file->Close();
}

void LArSlicer(const char* inputFile, const char* outputFile, int nmax = -1) {
    selectEvents(inputFile, outputFile, nmax);
}
