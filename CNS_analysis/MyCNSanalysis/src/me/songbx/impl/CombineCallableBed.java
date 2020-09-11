package me.songbx.impl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import me.songbx.model.FastaIndexEntry;

public class CombineCallableBed {

	/**
	 * The index are the names of transcripts, the values are the transcripts
	 */
	private HashMap<String, ArrayList<Integer>> callableHashMap = new HashMap<String, ArrayList<Integer>>();
	
	/**
	 * @param fileLocation
	 * @param chromoSomeRead
	 */
	
	public CombineCallableBed(String fastaFile, String bedListFile) {
		// read fasta file begin
		HashMap<String, FastaIndexEntry> entries = new ChromoSomeReadImplWithIndex(fastaFile).fastaIndexEntryImpl.getEntries();
		for ( String chr : entries.keySet() ) {
			callableHashMap.put(chr, new ArrayList<Integer>(Collections.nCopies(entries.get(chr).getLength(), 0)));			
		}
		// read fasta file end
		
		//read bed files begin
		File file = new File(bedListFile);
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            while ((tempString = reader.readLine()) != null) {
            	tempString = tempString.trim();
            	System.out.println(tempString);
            	// read bed file begin
            	GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(tempString));
            	BufferedReader br = new BufferedReader(new InputStreamReader(gzip));
            	
            	String tempString0 = null;
                while ((tempString0 = br.readLine()) != null) {
                	String[] arrOfStr = tempString0.split("\\s+");
                	if( arrOfStr[3].compareToIgnoreCase("CALLABLE") ==0 ){
                		int start = Integer.parseInt(arrOfStr[1]);
                		int end = Integer.parseInt(arrOfStr[2]);
                		for ( int i=start; i<end; ++i ) {
                			callableHashMap.get(arrOfStr[0]).set(i, callableHashMap.get(arrOfStr[0]).get(i) + 1);
                		}
                	}
                }
                br.close();
                // read bed file end
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {

                }
            }
        }
        //read bed files end
	}
	public HashMap<String, ArrayList<Integer>> getCallableHashMap() {
		return callableHashMap;
	}
	public void setCallableHashMap(HashMap<String, ArrayList<Integer>> callableHashMap) {
		this.callableHashMap = callableHashMap;
	}
	
	public static void main(String argv[]) {
		new CombineCallableBed("/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/reference.fa", "/media/bs674/panAndAssemblyfi/maize282Genotyping/callablePipeline/callable_list");
	}
}
