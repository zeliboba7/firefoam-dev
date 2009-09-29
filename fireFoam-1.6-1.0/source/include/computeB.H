
if (ignited) 
{
  b=b-(1.0/ignRampTime)*runTime.deltaT().value();
  b.max(0);
}
