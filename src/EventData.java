import java.io.Serializable;
import java.util.*;
public class EventData implements Serializable {

	public EventData()
	{
		pairData=new ArrayList<HadronPairData>();	
	}
	
	public boolean hasMC;
	public float Q2;
	public float W;
	public float x;
	public int beamPolarization;
	public ArrayList<HadronPairData> pairData;
	
}
