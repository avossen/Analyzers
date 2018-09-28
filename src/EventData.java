import java.io.Serializable;
import java.util.*;
public class EventData implements Serializable {

	public EventData()
	{
		pairData=new ArrayList<HadronPairData>();	
	}
	//the event and run number were not there before and probably break the schema
	//might be able to read old file if removed
	public int eventNr;
	public int runNr;
	public float torus;
	public float solenoid;
	public boolean hasMC;
	public float Q2;
	public float W;
	public float x;
	public float y;
	public float beamPolarization;
	public int beamHelicity;
	public ArrayList<HadronPairData> pairData;
	
}
