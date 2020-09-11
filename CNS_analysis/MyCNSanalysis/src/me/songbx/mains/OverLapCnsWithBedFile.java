package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.model.Cds;
import me.songbx.model.GeneSimple;
import me.songbx.model.Transcript;
import me.songbx.service.IfIntron;

public class OverLapCnsWithBedFile {
	public static void main(String[] args) {
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/sorghum/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.42.gff3");
		HashMap<String, ArrayList<Transcript>> transcriptHashMap = ifIntron.getTranscriptHashMap();
		HashMap<String, ArrayList<GeneSimple>> geneHashMap = ifIntron.getGeneHashMap();
		HashMap<String, HashSet<Integer>> allGeneticsBps = new HashMap<String, HashSet<Integer>>();
		HashMap<String, HashSet<Integer>> allCdsBps = new HashMap<String, HashSet<Integer>>();
		for( int i=1; i<=10; i++ ) {
			String chr = Integer.toString(i);
			allGeneticsBps.put(chr, new HashSet<Integer>());
			for( GeneSimple geneSimple : geneHashMap.get(chr) ) {
				for( int position=geneSimple.getStart(); position<=geneSimple.getEnd(); ++position ) {
					allGeneticsBps.get(chr).add(position);
				}
			}
			allCdsBps.put(chr, new HashSet<Integer>());
			for( Transcript transcript : transcriptHashMap.get(chr) ) {
				for( Cds cds :  transcript.getCdsHashSet()) {
					for( int position=cds.getStart(); position<=cds.getEnd(); ++position ) {
						allCdsBps.get(chr).add(position);
					}
				}
			}
		}
		// if you want to overlap with non-genic use allGeneticsBps,
		//if you want to overlap with non-coding then use allCdsBps
//		new OverLapCnsWithBedFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", allGeneticsBps, "/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed");
//		new OverLapCnsWithBedFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", allGeneticsBps, "/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed");
//		new OverLapCnsWithBedFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.wig", allGeneticsBps, "/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed");
//		new OverLapCnsWithBedFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.wig", allGeneticsBps, "/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed");
		
		new OverLapCnsWithBedFile("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.wig", allGeneticsBps, "/media/bs674/1_8t/AndCns/CNSexpression/transcript.exp_at.least.FPKM0.5_exon.bed");
		
    }
	
	public OverLapCnsWithBedFile(String wigFile, HashMap<String, HashSet<Integer>> excludeBps, String bedFile) {
		try {
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
			HashMap<String, HashSet<Integer>> bedRecords = new HashMap<String, HashSet<Integer>> ();
			HashMap<String, HashSet<Integer>> cnsRecords = new HashMap<String, HashSet<Integer>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				bedRecords.put(chr, new HashSet<Integer>());
				cnsRecords.put(chr, new HashSet<Integer>());
			}
			readBedFile(bedFile, bedRecords);
			readWigFile(wigFile, cnsRecords);
			double genomeSize = 0;
			double bedRegions = 0;
			double cnsRegions = 0;
			double cnsAndBedRegions = 0;
			
    		for( String chr : bedRecords.keySet() ) { // so here only chromosomes 1-10 were used
    			for( int i=0; i<maizeGenome.get(chr).length(); ++i ) {
    				if( excludeBps.get(chr).contains(i+1) ) {
    				}else {
    					++genomeSize;
						if( cnsRecords.get(chr).contains(i+1)  ) {
							if ( bedRecords.get(chr).contains(i+1)  ) {
		    					++cnsAndBedRegions;
		    				}else {
		    					++cnsRegions;
		    				}
						} else if ( bedRecords.get(chr).contains(i+1)  ) {
	    					++bedRegions;
	    				}
	    				
    				}
    			}
    		}
    		System.out.println("genomeSize=" + genomeSize);
    		System.out.println("cnsRegions=" + cnsRegions);
    		System.out.println("bedRegions=" + bedRegions);
    		System.out.println("cnsAndBedRegions=" + cnsAndBedRegions);
    		double fold = (cnsAndBedRegions/cnsRegions)/(bedRegions/genomeSize);
    		System.out.println("fold=" + fold);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public void readBedFile(String bedFile, HashMap<String, HashSet<Integer>> bedRecords) {
		try {
			File file = new File( bedFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				if( bedRecords.containsKey(arrOfStr[0]) ) {
					for( int i=Integer.parseInt(arrOfStr[1])+1; i<=Integer.parseInt(arrOfStr[2]); ++i ) {
						bedRecords.get(arrOfStr[0]).add(i);
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	public void readWigFile(String wigFile, HashMap<String, HashSet<Integer>> wigRecords) {
		try {
			File file = new File( wigFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String chr = "";
			Pattern p =Pattern.compile("chrom=(\\S+)");
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) == 't' ) {
					
				}else if( tempString.charAt(0) == 'v' ) {
					Matcher m = p.matcher(tempString);
					if( m.find() ) {
						chr = m.group(1);
					}
				}else {
					String[] arrOfStr = tempString.split("\\s+");
					if( wigRecords.containsKey(chr) ) {
						if( Integer.parseInt(arrOfStr[1]) > 0 ) {
							wigRecords.get(chr).add(Integer.parseInt(arrOfStr[0]));
						}
					}
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
//	
//	static boolean ifWithinBedFileRegion( int position, ArrayList<BedRecord> bedRecords ) {
//		int start = 0 ;
//		int end = bedRecords.size()-1;
//		int lastStart = start;
//		if( bedRecords.get(start).getStart() > position ){
//			return false;
//		}
//		if( bedRecords.get(end).getEnd() < position){
//			return false;
//		}
//		while(!((bedRecords.get(start).getStart() <= position) && (bedRecords.get(start+1).getEnd( )>= position))){
//			if((bedRecords.get(start).getEnd() < position)){
//				lastStart = start;
//				if(1 == (end - start)){
//					start = end;
//				}else{
//					start = (start+end)/2;
//				}
//			}else{
//				end = start;
//				start = lastStart;
//			}
//		}
//		if( bedRecords.get(start).getStart()<=position && position<=bedRecords.get(start).getEnd()) {
//			return true;
//		} else if ( bedRecords.get(start+1).getStart()<=position && position<=bedRecords.get(start+1).getEnd() ) {
//			return true;
//		}
//		return false;
//	}
}
