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
	public float y;
	public float beamPolarization;
	public int beamHelicity;
	public ArrayList<HadronPairData> pairData;
	
}
