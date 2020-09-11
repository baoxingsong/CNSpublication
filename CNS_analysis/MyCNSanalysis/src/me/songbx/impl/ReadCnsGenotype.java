package me.songbx.impl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.model.BedRecord;

public class ReadCnsGenotype {
	HashMap<String/*chr*/, HashMap<String/* cns id*/, HashMap<String/*accession id*/, Integer /*CNS length in this accession*/>>> cnsGenotype = 
			new HashMap<String, HashMap<String, HashMap<String, Integer>>>();
	HashMap<String/*CNS id*/, BedRecord > cnsRecords = new HashMap<String/*CNS id*/, BedRecord >();
	HashMap<String/*CNS id*/, Integer > cnsLength = new HashMap<String/*CNS id*/, Integer >();
	HashSet<String> allAccessions = new HashSet<String>();
	
	public HashMap<String, HashMap<String, HashMap<String, Integer>>> getCnsGenotype() {
		return cnsGenotype;
	}

	public void setCnsGenotype(HashMap<String, HashMap<String, HashMap<String, Integer>>> cnsGenotype) {
		this.cnsGenotype = cnsGenotype;
	}

	public HashMap<String, BedRecord> getCnsRecords() {
		return cnsRecords;
	}

	public void setCnsRecords(HashMap<String, BedRecord> cnsRecords) {
		this.cnsRecords = cnsRecords;
	}

	public HashMap<String, Integer> getCnsLength() {
		return cnsLength;
	}

	public void setCnsLength(HashMap<String, Integer> cnsLength) {
		this.cnsLength = cnsLength;
	}

	public ReadCnsGenotype( String fileLocation ) {
		try {
        	File file = new File(fileLocation);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = reader.readLine();
            String [] accessions = tempString.split("\t");
            for( int ei = 2; ei<accessions.length; ++ei ) {
            	allAccessions.add(accessions[ei]);
            }
            Pattern p = Pattern.compile("^(\\S+):(\\d+)\\-(\\d+)");
            while ((tempString = reader.readLine()) != null) {
            	String [] elements = tempString.split("\t");
            	assert( accessions.length == elements.length);
            	Matcher m = p.matcher(elements[0]);
            	if( m.find() ) {
            		String chr = m.group(1);
            		int start = Integer.parseInt(m.group(2));
            		int end = Integer.parseInt(m.group(3));
            		int length = Integer.parseInt(elements[1]);
            		BedRecord bedRecord = new BedRecord(start, end);
            		cnsRecords.put(elements[0], bedRecord);
            		cnsLength.put(elements[0], length);
            		
            		if( ! cnsGenotype.containsKey(chr) ) {
            			cnsGenotype.put(chr, new HashMap<String, HashMap<String, Integer>>());
            		}
            		cnsGenotype.get(chr).put(elements[0], new HashMap<String, Integer>());
            		for( int ei = 2; ei<accessions.length; ++ei ) {
            			cnsGenotype.get(chr).get(elements[0]).put(accessions[ei], Integer.parseInt(elements[ei]));
                	}
            	}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
	}
	
	public HashSet<String> getAllAccessions() {
		return allAccessions;
	}

	public void setAllAccessions(HashSet<String> allAccessions) {
		this.allAccessions = allAccessions;
	}

	public static void main( String argv[] ) {
		new ReadCnsGenotype("/home/bs674/cbsuem02/local/workdir/bs674/CNSgwas/all.CNS.genotype.txt");
		
	}
}
