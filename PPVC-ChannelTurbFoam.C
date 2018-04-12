/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    PPVC-ChannelTurbuFoam

Description
    Transient solver for incompressible flow using a pressure projection,
    rotational velocity-correction strategy.
    Following the paper "A robust and accurate outflow boundary condition for
    incompressible flow simulations on severely-truncated unbounded domains" by
    Dong et al.

    A constant pressure gradient to compensate for momentum loass due to shear
    stress is added to conserve mass flux.
\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "IFstream.H"
#include "OFstream.H"
//#include "LESfilter.H"
//#include "IOmanip.H" // for input/ouput format control
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "createGradP.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    scalar gamma0; // temporal discretization coefficient
    dimensionedScalar dt = runTime.deltaT();

    //U0.correctBoundaryConditions();

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName()<< endl;
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // Start with one first order backward step
        if(runTime.value() == dt.value())
        {
            Info<< "\n first time step: 1st order backward"
                    "temporal discretization applied\n" << endl;
            gamma0 = 1.0;
            Uhat   = U;
            Ustar  = U;
        } else {
            //Info<< "\nSecond-order backward temporal discretization\n"<< endl;
            gamma0 = 1.5;
            Uhat   = 2.*U - 0.5*U0;
            Ustar  = 2.*U - U0;
        }

        // Update in every time step
        dt = runTime.deltaT();  // update deltaT
        U0 = U;

        Info<< "Updating W^(n+1) velocity.... " << endl;
        W.correctBoundaryConditions();//update W =U^(n+1) at dirichlet boundary

        // explicit source terms
        nuEff = turbulence->nuEff();

        Q1 = nuEff * fvc::grad(fvc::div(U)());
        Q2 = fvc::grad(nuEff) & (fvc::grad(U));
        Q3 = fvc::div(nuEff*dev2(fvc::grad(U)().T()));

        // Pressure Poisson for new pressure----------------------------------//
        phi = (fvc::interpolate(U) & mesh.Sf());

        // Non-orthogonal pressure corrector loop
        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
            fvm::laplacian(p)
            ==
            fvc::div(
                Uhat/dt                                // temporal term
                - (fvc::div(phi, U) - fvc::div(phi)*U) // advection term
                + Q1 + Q2 + Q3        // explicit turbulent stress terms
                )
            );
            pEqn.setReference(pRefCell, pRefValue);

            if (nonOrth == nNonOrthCorr)
            {
                pEqn.solve(mesh.solver("pFinal"));
            }
            else
            {
                pEqn.solve(mesh.solver("p"));
            }
        }

        // update p at all boundaries
        Info << "Final Pressure Update....." << endl;
        p.correctBoundaryConditions();

        // Helmholz equation for velocity ------------------------------------//
        surfaceScalarField phiStar = (fvc::interpolate(Ustar) & mesh.Sf());
        fvVectorMatrix UEqn
            (
              fvm::Sp(gamma0/dt, U)
            - fvm::laplacian(nuEff, U)
            ==
            (
              Uhat/dt
            - (fvc::div(phiStar, Ustar) - fvc::div(phiStar)*Ustar)
            - fvc::grad(p)
            + Q3
            + flowDirection*gradP // additional source term
            )
        );
        UEqn.solve();

        #include "calculateGradP.H"

        #include "continuityErrs.H"

        turbulence->correct();

        runTime.write();

        #include "writeGradP.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
