package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.impl.LeadingSnps;
import me.songbx.service.IfIntron;
import me.songbx.service.IfNoCodingRna;

public class OverLapCnsWithEveryFunctionalThingForCNSClassification {
	
	
	class BedRecord implements Comparable<BedRecord>{
		private int start;
		private int end;
		private int length;
		private int tag;
		private String samRecords;
		// 0000 0001 CNS
		// 0000 0010 open chromatin
		// 0000 0100 TFBS
		public int getTag() {
			return tag;
		}
		public BedRecord(int start, int end, int length, String samRecords) {
			super();
			this.start = start;
			this.end = end;
			this.length=length;
			this.samRecords = samRecords;
			this.tag = 0;
		}
		public String getSamRecords() {
			return samRecords;
		}
		public void setSamRecords(String samRecords) {
			this.samRecords = samRecords;
		}
		public void setTag(int tag) {
			this.tag = tag;
		}
		public int getStart() {
			return start;
		}
		public void setStart(int start) {
			this.start = start;
		}
		public int getEnd() {
			return end;
		}
		public int getLength() {
			return length;
		}
		public void setLength(int length) {
			this.length = length;
		}
		public void setEnd(int end) {
			this.end = end;
		}
		public BedRecord(int start, int end) {
			super();
			this.start = start;
			this.end = end;
			this.tag = 0;
		}
		@Override
		public int compareTo(BedRecord arg0) {
			return this.getStart() - arg0.getStart();
		}
		public boolean overlapWith ( BedRecord bedRecord ) {
			if( (bedRecord.getStart() <=start && start <=bedRecord.getEnd()) ||
					(bedRecord.getStart() <=end && end <=bedRecord.getEnd()) ||
					(start <= bedRecord.getStart() && bedRecord.getStart() <= end) ||
					(start <= bedRecord.getEnd() && bedRecord.getEnd() <= end) 
					) {
				return true;
			}
			return false;
		}
	}
	
	public static void main(String[] args) {
		ArrayList<String> bamFiles = new ArrayList<String>();
		bamFiles.add("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5_uniq.sam");
//		bamFiles.add("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5_uniq.sam");
//		bamFiles.add("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5_uniq.sam");
//		bamFiles.add("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5_uniq.sam");
//		bamFiles.add("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5_uniq.sam");
//		
		
		new OverLapCnsWithEveryFunctionalThingForCNSClassification(bamFiles);
		ArrayList<String> bedFiles = new ArrayList<String>();
		bedFiles.add("/media/bs674/1_8t/AndCns/coreCNSBasedGwas/core_cns.bed");
		new OverLapCnsWithEveryFunctionalThingForCNSClassification(bedFiles, "bed");
		
    }
	
	public OverLapCnsWithEveryFunctionalThingForCNSClassification(ArrayList<String> samFiles) {
		// read eQTL result begin
		
		HashMap<String, HashSet<Integer>> leadingSnps = ((new LeadingSnps()).readSummary());
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		
		IfNoCodingRna ifNoCodingRna = new IfNoCodingRna("/media/bs674/1_8t/AndCns/noncodingRNA/jipb12708-sup-0003-supdata-s3.gtf");
		
		try {
			HashMap<String, ArrayList<BedRecord>> openChromatinRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tFBSBedRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> cnsRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> hi_cRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> h3k9acs = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tes = new HashMap<String, ArrayList<BedRecord>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				tFBSBedRecords.put(chr, new ArrayList<BedRecord>());
				openChromatinRecords.put(chr, new ArrayList<BedRecord>());
				hi_cRecords.put(chr, new ArrayList<BedRecord>());
				h3k9acs.put(chr, new ArrayList<BedRecord>());
				tes.put(chr, new ArrayList<BedRecord>());
				cnsRecords.put(chr, new ArrayList<BedRecord>());
			}
			
			
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_HiC_B73_leaf_merged_loops.txt", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/RNAPII.loop", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/h3k4me3.loop", hi_cRecords);
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/ATAC_seq_data/BEDfiles/Zm_anotation_rep1.intersect.rep2.acr.bed", openChromatinRecords);
			readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed", openChromatinRecords);
			readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398047_ATAC_B73_ear.filtered_ACR.bed", openChromatinRecords);
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/all_reproducible_peaks_summits_merged.bed/all_reproducible_peaks_summits_merged.bed", tFBSBedRecords);
			
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed", h3k9acs);
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed", h3k9acs);
			
			HashSet<String> teClasses = new HashSet<String>();
			teClasses.add("DHH");
			teClasses.add("RST");
			teClasses.add("DTT");
			teClasses.add("DTA");
			teClasses.add("DTC");
			readTeGff("/media/bs674/1_8t/AndCns/overlapCNSWithDIfferentTEsuperfamily/B73.structuralTEv2.fulllength.2018-09-19.gff3", tes, teClasses);
			

			for(String samFile : samFiles) {
				readBamFile(samFile, cnsRecords);
			}
//			System.out.println("bed files reading done");
			
			for( String chr : openChromatinRecords.keySet() ) {
				Collections.sort(cnsRecords.get(chr));
				Collections.sort(tFBSBedRecords.get(chr));
				Collections.sort(openChromatinRecords.get(chr));
				Collections.sort(hi_cRecords.get(chr));
				Collections.sort(h3k9acs.get(chr));
				Collections.sort(tes.get(chr));
			}
			
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/sorghumCnsClassification.sam");
//			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/allTheCnsClassification.sam");
			
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
						if ( openChromatin.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (cnsRecord.getTag() | 1) ); // open chromatin
							break;
						}
					}
					
					for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
						if ( tFBSBedRecord.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (cnsRecord.getTag() | 2) ); // TFBS
							break;
						}
					}
					
					if( leadingSnps.containsKey(chr) ) {
						for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
							if(leadingSnps.get(chr).contains(position)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 4) ); // leading eQTL
								break;
							}
						}
					}
					
					for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
						if(ifIntron.getElement(chr, position) == 1) { // not intergenetic
							cnsRecord.setTag( (cnsRecord.getTag() | 8) ); // intron
							break;
						}
					}
					
					for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
						if(ifIntron.getElement(chr, position) == 3) { // not intergenetic
							cnsRecord.setTag( (cnsRecord.getTag() | 16) ); // utr
							break;
						}
					}
					
					if ( ifNoCodingRna.isNocodingRna(  chr, cnsRecord.getStart(), cnsRecord.getEnd() )){
						cnsRecord.setTag( (cnsRecord.getTag() | 32) ); // non-coding rna
					}
					
					if( hi_cRecords.containsKey(chr) ) {
						for( BedRecord hi_cRecord : hi_cRecords.get(chr) ) {
							if(hi_cRecord.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 64) ); //loop
								break;
							}
						}
					}
					
					if( h3k9acs.containsKey(chr) ) {
						for( BedRecord h3k9ac : h3k9acs.get(chr) ) {
							if(h3k9ac.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 128) ); //enhancer
								break;
							}
						}
					}
					
					if( tes.containsKey(chr) ) {
						for( BedRecord te : tes.get(chr) ) {
							if(te.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 256) ); //enriched TEs
								break;
							}
						}
					}
				}
			}
			
			HashMap<Integer, Integer> cnsClassifies = new HashMap<Integer, Integer>();
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					int cnsClassify = cnsRecord.getTag();
					if( ( cnsClassify & 1) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 2) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 4) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 8) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 16) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 32) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 64) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 128) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 256) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					
					outPut.print(cnsRecord.getLength() + "\t" + cnsRecord.getSamRecords());
					
					if( cnsClassifies.containsKey(cnsRecord.getTag()) ) {
						cnsClassifies.put( cnsRecord.getTag(), cnsClassifies.get(cnsRecord.getTag())+1);
					}else {
						cnsClassifies.put( cnsRecord.getTag(), 1);
					}
				}
			}
			outPut.close();
			System.out.println("OC\tTFBS\teqtl\tintron\tUTR\tnon-coding-RNA\tHI-C\tenhancer\tTE\tcount");
			for( int cnsClassify : cnsClassifies.keySet()) {
				if( ( cnsClassify & 1) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 2) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 4) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 8) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 16) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 32) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 64) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 128) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 256) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				System.out.println(cnsClassifies.get(cnsClassify));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public OverLapCnsWithEveryFunctionalThingForCNSClassification(ArrayList<String> bedFiles, String bed /*this variable would not be used, just for different version of the same function*/) {
		// read eQTL result begin
		
		HashMap<String, HashSet<Integer>> leadingSnps = ((new LeadingSnps()).readSummary());
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		
		IfNoCodingRna ifNoCodingRna = new IfNoCodingRna("/media/bs674/1_8t/AndCns/noncodingRNA/jipb12708-sup-0003-supdata-s3.gtf");
		
		try {
			HashMap<String, ArrayList<BedRecord>> openChromatinRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tFBSBedRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> cnsRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> hi_cRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> h3k9acs = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tes = new HashMap<String, ArrayList<BedRecord>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				tFBSBedRecords.put(chr, new ArrayList<BedRecord>());
				openChromatinRecords.put(chr, new ArrayList<BedRecord>());
				hi_cRecords.put(chr, new ArrayList<BedRecord>());
				h3k9acs.put(chr, new ArrayList<BedRecord>());
				tes.put(chr, new ArrayList<BedRecord>());
				cnsRecords.put(chr, new ArrayList<BedRecord>());
			}
			
			
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_HiC_B73_leaf_merged_loops.txt", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/RNAPII.loop", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/h3k4me3.loop", hi_cRecords);
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/ATAC_seq_data/BEDfiles/Zm_anotation_rep1.intersect.rep2.acr.bed", openChromatinRecords);
			readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed", openChromatinRecords);
			readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398047_ATAC_B73_ear.filtered_ACR.bed", openChromatinRecords);
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/all_reproducible_peaks_summits_merged.bed/all_reproducible_peaks_summits_merged.bed", tFBSBedRecords);
			
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed", h3k9acs);
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed", h3k9acs);
			
			HashSet<String> teClasses = new HashSet<String>();
			teClasses.add("DHH");
			teClasses.add("RST");
			teClasses.add("DTT");
			teClasses.add("DTA");
			teClasses.add("DTC");
			readTeGff("/media/bs674/1_8t/AndCns/overlapCNSWithDIfferentTEsuperfamily/B73.structuralTEv2.fulllength.2018-09-19.gff3", tes, teClasses);
			

			for(String bedFile : bedFiles) {
				readBedFile(bedFile, cnsRecords);
			}
//			System.out.println("bed files reading done");
			
			for( String chr : openChromatinRecords.keySet() ) {
				Collections.sort(cnsRecords.get(chr));
				Collections.sort(tFBSBedRecords.get(chr));
				Collections.sort(openChromatinRecords.get(chr));
				Collections.sort(hi_cRecords.get(chr));
				Collections.sort(h3k9acs.get(chr));
				Collections.sort(tes.get(chr));
			}
			
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/coreCnsClassification.sam");
			
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
						if ( openChromatin.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (cnsRecord.getTag() | 1) ); // open chromatin
							break;
						}
					}
					
					for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
						if ( tFBSBedRecord.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (cnsRecord.getTag() | 2) ); // TFBS
							break;
						}
					}
					
					if( leadingSnps.containsKey(chr) ) {
						for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
							if(leadingSnps.get(chr).contains(position)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 4) ); // leading eQTL
								break;
							}
						}
					}
					
					for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
						if(ifIntron.getElement(chr, position) == 1) { // not intergenetic
							cnsRecord.setTag( (cnsRecord.getTag() | 8) ); // intron
							break;
						}
					}
					
					for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
						if(ifIntron.getElement(chr, position) == 3) { // not intergenetic
							cnsRecord.setTag( (cnsRecord.getTag() | 16) ); // utr
							break;
						}
					}
					
					if ( ifNoCodingRna.isNocodingRna(  chr, cnsRecord.getStart(), cnsRecord.getEnd() )){
						cnsRecord.setTag( (cnsRecord.getTag() | 32) ); // non-coding rna
					}
					
					if( hi_cRecords.containsKey(chr) ) {
						for( BedRecord hi_cRecord : hi_cRecords.get(chr) ) {
							if(hi_cRecord.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 64) ); //loop
								break;
							}
						}
					}
					
					if( h3k9acs.containsKey(chr) ) {
						for( BedRecord h3k9ac : h3k9acs.get(chr) ) {
							if(h3k9ac.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 128) ); //enhancer
								break;
							}
						}
					}
					
					if( tes.containsKey(chr) ) {
						for( BedRecord te : tes.get(chr) ) {
							if(te.overlapWith(cnsRecord)) {
								cnsRecord.setTag( (cnsRecord.getTag() | 256) ); //enriched TEs
								break;
							}
						}
					}
				}
			}
			
			HashMap<Integer, Integer> cnsClassifies = new HashMap<Integer, Integer>();
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					int cnsClassify = cnsRecord.getTag();
					if( ( cnsClassify & 1) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 2) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 4) > 0){				
						outPut.print("1\t");
					}else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 8) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 16) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 32) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 64) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 128) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					if( ( cnsClassify & 256) > 0){				
						outPut.print("1\t");
					} else {
						outPut.print("0\t");
					}
					
					outPut.print(cnsRecord.getLength() + "\t" + cnsRecord.getSamRecords());
					
					if( cnsClassifies.containsKey(cnsRecord.getTag()) ) {
						cnsClassifies.put( cnsRecord.getTag(), cnsClassifies.get(cnsRecord.getTag())+1);
					}else {
						cnsClassifies.put( cnsRecord.getTag(), 1);
					}
				}
			}
			outPut.close();
			System.out.println("OC\tTFBS\teqtl\tintron\tUTR\tnon-coding-RNA\tHI-C\tenhancer\tTE\tcount");
			for( int cnsClassify : cnsClassifies.keySet()) {
				if( ( cnsClassify & 1) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 2) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 4) > 0){				
					System.out.print("1\t");
				}else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 8) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 16) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 32) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 64) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 128) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				if( ( cnsClassify & 256) > 0){				
					System.out.print("1\t");
				} else {
					System.out.print("0\t");
				}
				System.out.println(cnsClassifies.get(cnsClassify));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void readBedFile(String bedFile, HashMap<String, ArrayList<BedRecord>> bedRecords) {
		try {
			File file = new File( bedFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				if( bedRecords.containsKey(arrOfStr[0]) ) {
					int start = Integer.parseInt(arrOfStr[1])+1;
					int end = Integer.parseInt(arrOfStr[2]);
					BedRecord bedRecord = new BedRecord(start, end);
					bedRecords.get(arrOfStr[0]).add(bedRecord);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readTeGff(String gffFile, HashMap<String, ArrayList<BedRecord>> bedRecords, HashSet<String> teClasses) {
		try {
			File file = new File( gffFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			Pattern p = Pattern.compile("(\\S+)\\s+\\S+\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+.*ID=((\\w{3})\\d+)");
			while ((tempString = reader.readLine()) != null) {
				Matcher m = p.matcher(tempString);
				if( m.find() && bedRecords.containsKey(m.group(1)) && teClasses.contains(m.group(5)) ) {
					int start = Integer.parseInt(m.group(2));
					int end = Integer.parseInt(m.group(3));
					BedRecord bedRecord = new BedRecord(start, end);
					bedRecords.get(m.group(1)).add(bedRecord);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readBamFile(String bamFile, HashMap<String, ArrayList<BedRecord>> bedRecords) {
		try {
			SamReader reader0 = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
			SAMRecordIterator it = reader0.iterator();
			while (it.hasNext()) {
	            final SAMRecord samRecord = it.next();
	            int start = samRecord.getAlignmentStart();
	            int end = samRecord.getAlignmentEnd();
	            String chr = samRecord.getContig();
	            if( bedRecords.containsKey(chr) ) {
	            	int matchLength = 0;
	            	Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
    				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
    				while (matcher.find()) {
    					int length = Integer.parseInt( matcher.group(1));
    					char cigerLetter =  matcher.group(2).charAt(0);
    					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
    						matchLength += length;
    					}
    				}
	            	BedRecord bedRecord = new BedRecord(start, end, matchLength, samRecord.getSAMString());
	            	bedRecords.get(chr).add(bedRecord);
	            }
			}
			reader0.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void readDnaLoopFile(String bedFile, HashMap<String, ArrayList<BedRecord>> bedRecords) {
		try {
			File file = new File( bedFile );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				String[] arrOfStr = tempString.split("\\s+");
				if( bedRecords.containsKey(arrOfStr[0]) ) {
					int start = Integer.parseInt(arrOfStr[1])+1;
					int end = Integer.parseInt(arrOfStr[2]);
					BedRecord bedRecord = new BedRecord(start, end);
					bedRecords.get(arrOfStr[0]).add(bedRecord);
				}
				if( bedRecords.containsKey(arrOfStr[3]) ) {
					int start = Integer.parseInt(arrOfStr[4])+1;
					int end = Integer.parseInt(arrOfStr[5]);
					BedRecord bedRecord = new BedRecord(start, end);
					bedRecords.get(arrOfStr[0]).add(bedRecord);
				}
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
