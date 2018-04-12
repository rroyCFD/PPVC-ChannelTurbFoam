Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U0
(
    IOobject
    (
        "U0",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);   

volVectorField W
(
    IOobject
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);       
    //bool U0Exist = true;
    
    // volVectorField W
    // (
    //     IOobject
    //     (
    //         "W",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );    

//---------------------------------------------------------------------------------//    
//     IOobject U0Header
//     (
//         "U0",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ
//     );
//     
//     volVectorField* U0Pointer;
//     bool U0Exist = false;
// 
//     if (U0Header.headerOk())
//     {
//         U0Exist = true;
//         Info<< "Reading field U0\n" << endl;
//         
//         U0Pointer = new volVectorField
//         (
//             IOobject
//             (
//                 "U0",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh
//         );
//     } else {
//         Info<< "U0 unavaiable! Initializing U0 = U \n" << endl;
//         U0Pointer = new volVectorField
//         (
//             IOobject
//             (
//                 "U0",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             //U 
//             //mesh,
//             //dimensionedVector("U0",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero),
//             //U.boundaryField().types()
//         );  
//     }    
//     
//     volVectorField U0 = *U0Pointer;
//---------------------------------------------------------------------------------//
    
Info<< "Creating field Ustar \n" << endl;
volVectorField Ustar
(
    IOobject
    (
        "Ustar",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("Ustar",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
);

Info<< "Creating field Uhat \n" << endl;
volVectorField Uhat
(
    IOobject
    (
        "Uhat",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("Uhat",dimensionSet(0, 1, -1, 0, 0, 0, 0),vector::zero)
);

    
#include "createPhi.H"   
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


//------------------------experiemntal branch --------------------------------//
volScalarField nuEff
(
    IOobject
    (
        "nuEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    //turbulence->nuEff() // can't find in lookupObjectRegistry, in BC implementation
    mesh,
    dimensionedScalar("nuEff",dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
);
nuEff = turbulence->nuEff();

volVectorField Q1
(
    IOobject
    (
        "Q1",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("Q1",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
);

volVectorField Q2
(
    IOobject
    (
        "Q2",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("Q2",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
);


volVectorField Q3
(
    IOobject
    (
        "Q3",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedVector("Q3",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
);
//----------------------------------------------------------------------------//


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
        
    
const dictionary& ProjectionSolverDict = mesh.solutionDict().subDict("ProjectionSolver");
//const int nCorr = ProjectionSolverDict.lookupOrDefault<int>("nCorrectors", 1);
const int nNonOrthCorr = ProjectionSolverDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 1);

const bool pRefOn = ProjectionSolverDict.lookupOrDefault("pRefOn", false);    
label  pRefCell =  ProjectionSolverDict.lookupOrDefault<scalar>("pRefCell", 0); // label pRefCell = 0;
scalar pRefValue = ProjectionSolverDict.lookupOrDefault<scalar>("pRefValue", 0); // scalar pRefValue = 0.0;

// definition: https://github.com/OpenFOAM/OpenFOAM-2.2.x/blob/master/src/finiteVolume/cfdTools/general/findRefCell/findRefCell.C
// search for RefCell and RefValue within. So initializing with an arbirary value also works!
setRefCell(p, ProjectionSolverDict, pRefCell, pRefValue); 

Info << "pRefOn: " << pRefOn << endl;
if(pRefOn) {
    Info <<"pRefCell: "<< pRefCell << "  pRefValue: "<< pRefValue <<endl;
}
          
// // LES Filter 
// autoPtr<LESfilter> filterPtr(LESfilter::New(U.mesh(), transportProperties));
// LESfilter& filter(filterPtr());

