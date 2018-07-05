import java.io.File;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.io.*;
import java.nio.file.*;
import java.io.Serializable;
import java.util.*;


public class PairReader {

	public static void main(String[] args) 
	{
		// TODO Auto-generated method stub
		AsymData m_asymData = null;
		File folder = new File(args[0]);
		File[] listOfFiles = folder.listFiles();
		for (int iF = 0; iF < listOfFiles.length; iF++) {
			if (listOfFiles[iF].isFile()) {
				System.out.println("File " + listOfFiles[iF].getName());
				PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:*.{srn}");

				Path filename = Paths.get(listOfFiles[iF].getName());
				if (matcher.matches(filename)) 
				{
	        // Deserialization
					try
					{   
	            // Reading the object from a file
						FileInputStream file = new FileInputStream(filename.toString());
						ObjectInputStream in = new ObjectInputStream(file);
	             
	            // Method for deserialization of object
						m_asymData = (AsymData)in.readObject();
						System.out.println("got "+ m_asymData.eventData.size() + " hadron pairs");
						for(int i =0;i<m_asymData.eventData.size();i++)
						{
							for(int j=0;j<m_asymData.eventData.get(i).pairData.size();j++)
							{
								System.out.println("looking at hadron pair with z: "+ m_asymData.eventData.get(i).pairData.get(j).z);			
							}
						}
						in.close();
						file.close();
	             
						System.out.println("Object has been deserialized "); 
	        }	         
	        catch(IOException ex)
	        {
	            System.out.println("IOException is caught");
	        }
	        catch(ClassNotFoundException ex)
	        {
	            System.out.println("ClassNotFoundException is caught");
	        }
	 
	    }
				}
			}
	}	
	}
