// Inside the time loop, after solving:
if (runTime.timeIndex() % checkInterval == 0)
{
    scalar maxDiff = max(mag(alpha1 - alpha1Old)).value();
    
    Info << "Interface movement metric: " << maxDiff << endl;

    // check max value of alphaPhi10
    Info << "alphaPhi10 mean: " << average(alphaPhi10) << endl;
    Info << "alphaPhi10 max: " << max(alphaPhi10) << endl;
    
    if (maxDiff < stopTolerance)
    {
        Info << "Steady state reached at t=" 
             << runTime.timeName() << endl;
        runTime.writeAndEnd();  // writes final fields and stops cleanly
    }
    
    alpha1Old = alpha1;  // update reference field

}

