{
    // Correct driving force for a constant mass flow rate
    volScalarField rAU(1.0/UEqn.A());

    // Extract the mean velocity in the flow direction
    dimensionedScalar magUbarStar =
        (flowDirection & U)().weightedAverage(mesh.V());

    // Calculate the pressure gradient increment needed to
    // adjust the average flow-rate to the correct value
    dimensionedScalar gradPplus =
        (magUbar - magUbarStar)/rAU.weightedAverage(mesh.V());

    U += flowDirection*rAU*gradPplus;
    U.correctBoundaryConditions();

    gradP += gradPplus;

    Info<< "Uncorrected Ubar = " << magUbarStar.value() << tab
        << "pressure gradient = " << gradP.value() << endl;
}
