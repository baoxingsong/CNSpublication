package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.service.IfIntron;

public class OverLapCnsWithTFBSOpenChromatinForCNSClassification {
	
	
	class BedRecord implements Comparable<BedRecord>{
		private int start;
		private int end;
		private int tag;
		private String samRecords;
		// 0000 0001 CNS
		// 0000 0010 open chromatin
		// 0000 0100 TFBS
		public int getTag() {
			return tag;
		}
		public BedRecord(int start, int end, String samRecords) {
			super();
			this.start = start;
			this.end = end;
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
		bamFiles.add("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5.bam");

		//		bamFiles.add("/media/bs674/1_8t/AndCns/Setaria_italica/result/5.bam");
		new OverLapCnsWithTFBSOpenChromatinForCNSClassification(bamFiles);
		
    }
	
	public OverLapCnsWithTFBSOpenChromatinForCNSClassification(ArrayList<String> samFiles) {
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		try {
			HashMap<String, ArrayList<BedRecord>> openChromatinRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tFBSBedRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> enhancers = new HashMap<String, ArrayList<BedRecord>> ();
			
			HashMap<String, ArrayList<BedRecord>> cnsRecords = new HashMap<String, ArrayList<BedRecord>> ();

			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				tFBSBedRecords.put(chr, new ArrayList<BedRecord>());
				openChromatinRecords.put(chr, new ArrayList<BedRecord>());
				cnsRecords.put(chr, new ArrayList<BedRecord>());
				enhancers.put(chr, new ArrayList<BedRecord>());
			}
			
			ArrayList<String> openChromatinFiles = new ArrayList<String>();
			openChromatinFiles.add("/media/bs674/2t/openChromotainAndTFBS/ATAC_seq_data/BEDfiles/Zm_anotation_rep1.intersect.rep2.acr.bed");
			openChromatinFiles.add("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed");
			openChromatinFiles.add("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398047_ATAC_B73_ear.filtered_ACR.bed");
			for( String openChromatinFile : openChromatinFiles ) {
				readBedFile(openChromatinFile, openChromatinRecords);
			}
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/all_reproducible_peaks_summits_merged.bed/all_reproducible_peaks_summits_merged.bed", tFBSBedRecords);
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed", enhancers);
			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed", enhancers);
			
			PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEnhancerCns.sam");
			PrintWriter outPut1 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBSOrEnhancer.sam");
			PrintWriter outPut2 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsIntronUtrCns.sam");
			PrintWriter outPut3 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.sam");
			PrintWriter outPut4 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam");
			
			PrintWriter outPut5 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntron_openchromatin.sam");
			PrintWriter outPut6 = new PrintWriter("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_noneopenchromatin.sam");
			
			for(String samFile : samFiles) {
				readBamFile(samFile, cnsRecords);
			}
			System.out.println("bed files reading done");
			
			for( String chr : openChromatinRecords.keySet() ) {
				Collections.sort(cnsRecords.get(chr));
				Collections.sort(tFBSBedRecords.get(chr));
				Collections.sort(openChromatinRecords.get(chr));
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					cnsRecord.setTag( (short)1 );
				}
			}
			
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
						if ( openChromatin.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)2) );
							break;
						}
					}
					
					for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
						if ( tFBSBedRecord.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)4) );
							break;
						}
					}
					
					for( BedRecord enhancer : enhancers.get(chr) ) {
						if ( enhancer.overlapWith(cnsRecord) ){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)8) );
							break;
						}
					}
					
					for( int position = cnsRecord.getStart(); position<=cnsRecord.getEnd(); ++position ) {
						if(ifIntron.getElement(chr, position) != 0) { // not intergenetic
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)16) );
							break;
						}
					}
				}
			}
			HashMap<Integer, Integer> cnsClassifies = new HashMap<Integer, Integer>();
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					if( cnsRecord.getTag() == 1 || cnsRecord.getTag() == 17 ) { // 17 did not overlap with enhancer, open chromatin and TFBS, but overlapped with genetic region
						outPut.print(cnsRecord.getSamRecords());
					}
					if( ((cnsRecord.getTag() & (short)16) == 0) && ((cnsRecord.getTag() & (short)2) == 0) && ((cnsRecord.getTag() & (short)4) == 0) ) {  
						outPut2.print(cnsRecord.getSamRecords());
					}
					
					if( ((cnsRecord.getTag() & (short)16) == 0) && ( ((cnsRecord.getTag() & (short)2) != 0) || ((cnsRecord.getTag() & (short)4) != 0) ) ) {  
						outPut4.print(cnsRecord.getSamRecords());
					}
					
					if( ((cnsRecord.getTag() & (short)16) == 0) && ( ((cnsRecord.getTag() & (short)2) != 0) || ((cnsRecord.getTag() & (short)4) != 0)
							|| ((cnsRecord.getTag() & (short)8) != 0) ) ) {  
						outPut1.print(cnsRecord.getSamRecords());
					}
					
					if( (cnsRecord.getTag() & (short)16) == 0 ) {  // did not overlapped with genetic region
						outPut3.print(cnsRecord.getSamRecords());
					}
					if( cnsClassifies.containsKey((int)cnsRecord.getTag()) ) {
						cnsClassifies.put((int) cnsRecord.getTag(), cnsClassifies.get((int)cnsRecord.getTag())+1);
					}else {
						cnsClassifies.put((int) cnsRecord.getTag(), 1);
					}
					if( (cnsRecord.getTag() & (short)2) == 0 && (cnsRecord.getTag() & (short)16) == 0 ) {  // did not overlapped with open chromatin
						outPut6.print(cnsRecord.getSamRecords());
					} else if ((cnsRecord.getTag() & (short)16) == 0) {
						outPut5.print(cnsRecord.getSamRecords());
					}
				}
			}
			outPut.close();
			outPut1.close();
			outPut2.close();
			outPut3.close();
			outPut4.close();
			outPut5.close();
			outPut6.close();
			for( int cnsClassify : cnsClassifies.keySet()) {
				System.out.println(cnsClassify + "\t" + cnsClassifies.get(cnsClassify));
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
	            	BedRecord bedRecord = new BedRecord(start, end, samRecord.getSAMString());
	            	bedRecords.get(chr).add(bedRecord);
	            }
			}
			reader0.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
