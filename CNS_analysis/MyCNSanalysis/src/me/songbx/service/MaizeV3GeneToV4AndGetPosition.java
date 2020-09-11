package me.songbx.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

public class MaizeV3GeneToV4AndGetPosition {
	HashMap<String, Integer> v3Counts = new HashMap<String, Integer>();
	HashMap<String, String> v3ToV4 = new HashMap<String, String>();	
	public MaizeV3GeneToV4AndGetPosition( String inputFile ) {
		try {
        	File file = new File(inputFile);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            while ((tempString = reader.readLine()) != null) {
            	String [] elements = tempString.split("\t");
            	
            	if( v3Counts.containsKey(elements[0]) ){
            		v3Counts.put(elements[0], v3Counts.get(elements[0])+1);
	            } else {
	            	v3Counts.put(elements[0], 1);
	            }
            	v3ToV4.put(elements[0], elements[1]);
            }
            reader.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
	}
	public HashMap<String, Integer> getV3Counts() {
		return v3Counts;
	}
	public void setV3Counts(HashMap<String, Integer> v3Counts) {
		this.v3Counts = v3Counts;
	}
	public HashMap<String, String> getV3ToV4() {
		return v3ToV4;
	}
	public void setV3ToV4(HashMap<String, String> v3ToV4) {
		this.v3ToV4 = v3ToV4;
	}
	public static void main( String argv[] ) {
		new MaizeV3GeneToV4AndGetPosition("/media/bs674/1_8t/AndCns/subfunctionlization/sorghum/correlationCNSWithTissueSpecificExpression/v3_to_v4");
	}
}
