import java.util.ArrayList;

public enum Binning {
MBinning(0),
ZBinning(1),
XBinning(2);

	
	protected int binType;
	protected ArrayList<Double> bins;
	public static final int numKinBins=3;
	
	Binning(int iBin)
	{
		bins=new ArrayList<Double >();
		binType=iBin;
		//m binning
		if(iBin==0)
		{
	//		bins.add(0.4);
			//bins.add(0.8);
			bins.add(1000.0);
		
		}
		//z binning
		if(iBin==1)
		{
	//		bins.add(0.3);
			//bins.add(0.6);
			bins.add(1.2);
		}
		//x binning
		if(iBin==2)
		{
	//		bins.add(0.05);
			//bins.add(0.2);
			bins.add(1.2);	
		}
	}
	
	public int getNumBins()
	{
		return bins.size();
	}
	//unnecessary since enum already has name() function
	public String getBinningName()
	{
		if(binType==0)
			return "M";
		if(binType==1)
			return "Z";
		if(binType==2)
			return "X";
					
		return "none";	
	}
	public String getBinningName(int ibin)
	{
		if(ibin==0)
			return "M";
		if(ibin==1)
			return "Z";
		if(ibin==2)
			return "X";
					
		return "none";	
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
		
		//System.out.println("getting bin for binType " + binType + " value: " + value);
	  int coo1=-1;

	  for(int i=0;i<bins.size();i++)
	    {
		//  System.out.println("compare to " + bins.get(i));
	      if(value<=bins.get(i))
	      {
	    	  	coo1=i;
	    	  	break;
	      }
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
      {
	coo1=i;
	break;
      }
    }
  /*  if(coo1<0)
    {
        cout <<"wrong coo: val: " << value <<endl;
	}*/
  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
  return coo1;
}
	
}
