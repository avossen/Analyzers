import java.util.ArrayList;

public enum Binning {
MBinning(0),
ZBinning(1),
XBinning(2),none(3);

	
	protected int binType;
	protected ArrayList<Double> bins;
	public final int numKinBins=3;
	
	Binning(int iBin)
	{
		binType=iBin;
		//m binning
		if(iBin==0)
		{
			bins.add(0.4);
			bins.add(0.8);
			bins.add(1000);
		
		}
		//z binning
		if(iBin==1)
		{
			bins.add(0.3);
			bins.add(0.6);
			bins.add(1.2);
		}
		//x binning
		if(iBin==2)
		{
			bins.add(0.05);
			bins.add(0.2);
			bins.add(1.2);	
		}
	}
	
	public int getBin(double value1,double value2, double value3)
	{
		double value=0;
		if(binType==0)
		{
			value=value1;
		}
		if(binType==1)
			value=value2;
		if(binType==2)
			value=value3;
		
	  int coo1=-1;

	  for(int i=0;i<bins.size();i++)
	    {
	      if(value<=bins.get(i))
	    
	    	  coo1=i;
	      break;
		}
	    
	  /*  if(coo1<0)
	    {
	        cout <<"wrong coo: val: " << value <<endl;
		}*/
	  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
	  return coo1;
	}

public int getBinType()
{
return binType;	
}
//for the phi bins(external)
public int getBin(ArrayList<Double> b1, double value)
{
  int coo1=-1;

  for(int i=0;i<b1.size();i++)
    {
      if(value<=b1.get(i))
    
	coo1=i;
	break;
	
    }
  /*  if(coo1<0)
    {
        cout <<"wrong coo: val: " << value <<endl;
	}*/
  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
  return coo1;
}
	
}
