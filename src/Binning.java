import java.util.ArrayList;

public enum Binning {
MBinning(0),
ZBinning(1),
XBinning(2),none(3);

	protected ArrayList<Double> bins;
	public final int numKinBins=3;
	
	Binning(int iBin)
	{
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
	
	
	
}
