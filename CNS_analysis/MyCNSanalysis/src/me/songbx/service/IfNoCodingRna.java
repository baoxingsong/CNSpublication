package me.songbx.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.model.BedRecord;

public class IfNoCodingRna {
	private HashMap<String, ArrayList<BedRecord>> nonCodingGeneHashMap;
	public IfNoCodingRna( String fileLocation ) {
		nonCodingGeneHashMap = new HashMap<String, ArrayList<BedRecord>>();
		HashSet<String> noCodingKeyWords = new HashSet<String>();
		noCodingKeyWords.add("lincRNA");
//		noCodingKeyWords.add("lincRNA_gene");
		noCodingKeyWords.add("miRNA");
//		noCodingKeyWords.add("miRNA_gene");
		noCodingKeyWords.add("tRNA_gene");
		
		File file = new File(fileLocation);
		BufferedReader reader = null;
        try {
        	reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            Pattern p = Pattern.compile("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+");
			while ((tempString = reader.readLine()) != null) {
            	Matcher m=p.matcher(tempString);
				if(tempString.startsWith("#")){
					
				}else{
					if(m.find()){
	                	if( noCodingKeyWords.contains(m.group(3)) || (m.group(2).compareTo("Cufflinks")==0 && m.group(3).compareTo("exon")==0 ) ){
	                		int start = Integer.parseInt(m.group(4));
		                	int end = Integer.parseInt(m.group(5));
		                	if(start>end){
		                		int temp=start;
		                		start = end;
		                		end = temp;
		                	}
		                	String chrId = m.group(1);
	                		BedRecord bedRecord = new BedRecord(start, end);
	                		if( ! nonCodingGeneHashMap.containsKey(chrId) ) {
	                			nonCodingGeneHashMap.put(chrId, new ArrayList<BedRecord>());
	                		}
	                		nonCodingGeneHashMap.get(chrId).add(bedRecord);
	                	}
	                }
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (reader != null) {
                try {
                    reader.close();
                } catch (IOException e1) {
                	e1.getStackTrace();
                }
            }
        }
        for ( String chr : nonCodingGeneHashMap.keySet() ) {
			Collections.sort(nonCodingGeneHashMap.get(chr));
		}
	}
	
	public boolean isNocodingRna ( String chr, int position1, int position2 ) {
		if( this.nonCodingGeneHashMap.containsKey(chr) ) {
			for( BedRecord nonCoding : nonCodingGeneHashMap.get(chr) ) {
				if( (nonCoding.getStart() <=position1 && position1 <=nonCoding.getEnd()) ||
						(nonCoding.getStart() <=position2 && position2 <=nonCoding.getEnd()) ||
						(position1 <= nonCoding.getStart() && nonCoding.getStart() <= position2) ||
						(position1 <= nonCoding.getEnd() && nonCoding.getEnd() <= position2) 
						) {
					return true;
				}
			}
		}
		return false;
	}
	public static void main( String argv[]) {
		IfNoCodingRna ifNoCodingRna = new IfNoCodingRna("/media/bs674/1_8t/AndCns/noncodingRNA/jipb12708-sup-0003-supdata-s3.gtf");
		if( ifNoCodingRna.isNocodingRna("1", 572435, 572438) ) {
			System.out.println("yes");
		}else {
			System.out.println("no");
		}
	}
}
