import java.io.Serializable;
import java.util.*;

public class AsymData implements Serializable {

	//has associated MC data
		public boolean hasMC;
		
	//
		public ArrayList<EventData> eventData;
		public ArrayList<EventData> eventDataMC;
		public AsymData()
		{
			eventData=new ArrayList<EventData>();
			eventDataMC=new ArrayList<EventData>();
			
		}
}
