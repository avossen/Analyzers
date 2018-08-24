import java.io.Serializable;
import java.util.*;

public class HadronPairData implements Serializable {

	public float z;
	public float xF;
	public float weight;
	public float phiR;
	public float phiH;
	public float M;
	public float theta;
	public boolean hasMC;
	public float pTLab;
	public float pTBreit;
	public float theta1;
	public float theta2;
	public float phi1;
	public float phi2;
	public float mom1;
	public float mom2;
	
	HadronPairData matchingMCPair;
	
	
	public HadronPairData()
	{
		
	}
	
	
}
