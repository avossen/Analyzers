import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;
import NovelFitters.MyParticle;

public class HadronPair {

	protected double theta;
	protected double phi_h;
	
	protected double z;
	protected double xF;
	
	protected Vector3 Ph;
	protected double pT;
	
	public double getTheta()
	{
		return theta;	
	}
	public double getPhiH()
	{
		return phi_h;
	}
	public double getZ()
	{
		return z;
	}
	
	public double getXf()
	{
		return xF;
	}
	
	public Vector3 geetPh()
	{
		return Ph;
	}
	public double getPt()
	{
		return pT;
	}
	
	//takes virtual photon  and Was parameter
	HadronPair(MyParticle h1, MyParticle h2, LorentzVector lv_q,double W)
	{
		compute();	
	
	}
	
	public void compute()
	{
		
	
	}
	
	
}
