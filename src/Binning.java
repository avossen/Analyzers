import java.util.ArrayList;

public enum Binning {
MBinning(0),
ZBinning(1),
XBinning(2),
RunBinning(3);
	
	protected int binType;
	protected ArrayList<Double> bins;
	public static final int numKinBins=4;
	
	Binning(int iBin)
	{
		bins=new ArrayList<Double >();
		binType=iBin;
		//m binning
		if(iBin==0)
		{
			bins.add(0.4);
			bins.add(0.8);
			bins.add(1000.0);
		
		}
		//z binning
		if(iBin==1)
		{
			bins.add(0.4);
			bins.add(0.6);
			bins.add(1.2);
		}
		//x binning
		if(iBin==2)
		{
		
			bins.add(0.2);
			bins.add(0.4);
			bins.add(1.2);	
		}
		//runNumber
		if(iBin==3)
		{
			//outbending
			for(int i=3910;i<3978;i++)
			{
				bins.add((double)i);
			}
			
			
			bins.add(4013.0);
			bins.add(4014.0);
			bins.add(4015.0);
			bins.add(4016.0);
			bins.add(4017.0);
			bins.add(4018.0);
			bins.add(4020.0);
			bins.add(4021.0);
			bins.add(4022.0);
			bins.add(4025.0);
			bins.add(4026.0);
			bins.add(4027.0);
			bins.add(4028.0);
			bins.add(4030.0);
			bins.add(4032.0);
			bins.add(4033.0);
			bins.add(4037.0);
			bins.add(4038.0);
			bins.add(4039.0);
			bins.add(4041.0);
			bins.add(4044.0);
			bins.add(4050.0);
			bins.add(4053.0);
			bins.add(4060.0);
			bins.add(4061.0);
			bins.add(4067.0);
			bins.add(4068.0);
			bins.add(4069.0);
			bins.add(4070.0);
			bins.add(4071.0);
			bins.add(4073.0);
			bins.add(4074.0);
			bins.add(4075.0);
			bins.add(4078.0);
			bins.add(4301.0);
			bins.add(4302.0);
			bins.add(4303.0);
			bins.add(4304.0);
			bins.add(4305.0);
			bins.add(4306.0);
			bins.add(4307.0);
			bins.add(4308.0);
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
		if(binType==3)
			return "run";	
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
		if(ibin==3)
			return "run";
		return "none";	
	}
	
	
	public int getBin(double value1,double value2, double value3, double value4)
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
		if(binType==3)
			value=value4;
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
public static int getBin(ArrayList<Double> b1, double value)
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
