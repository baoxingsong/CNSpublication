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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.impl.ChromoSomeReadImpl;

public class OverLapCnsWithTFBSOpenChromatinEQTLForKmer {
	class BedRecord implements Comparable<BedRecord>{
		private int start;
		private int end;
		private short tag;
		// 0000 0001 CNS
		// 0000 0010 open chromatin
		// 0000 0100 TFBS
		String querySequence;
		public String getQuerySequence() {
			return querySequence;
		}
		public void setQuerySequence(String querySequence) {
			this.querySequence = querySequence;
		}
		public short getTag() {
			return tag;
		}
		public void setTag(short tag) {
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
		new OverLapCnsWithTFBSOpenChromatinEQTLForKmer("/media/bs674/1_8t/AndCns/CNSprofiles/sorghum_nonIntronUtrCns.sam");
    }
	public OverLapCnsWithTFBSOpenChromatinEQTLForKmer(String samFile) {
		try {
			int seq_length = 50;
			HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
			HashMap<String, String> maizeGenome2 = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46_gene.fa").getChromoSomeHashMap();
			
			HashMap<String, ArrayList<BedRecord>> openChromatinRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> tFBSBedRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> cnsRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> hi_cRecords = new HashMap<String, ArrayList<BedRecord>> ();
			HashMap<String, ArrayList<BedRecord>> h3k9acs = new HashMap<String, ArrayList<BedRecord>> ();
			
			for( int i=1; i<=10; ++i ) {
				String chr = Integer.toString(i);
				tFBSBedRecords.put(chr, new ArrayList<BedRecord>());
				openChromatinRecords.put(chr, new ArrayList<BedRecord>());
				hi_cRecords.put(chr, new ArrayList<BedRecord>());
				h3k9acs.put(chr, new ArrayList<BedRecord>());
				cnsRecords.put(chr, new ArrayList<BedRecord>());
			}
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/ATAC_seq_data/BEDfiles/Zm_anotation_rep1.intersect.rep2.acr.bed", openChromatinRecords);
			//readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398046_ATAC_B73_leaf.filtered_ACR.bed", openChromatinRecords);
			//readBedFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_RAW/GSM3398047_ATAC_B73_ear.filtered_ACR.bed", openChromatinRecords);
			
			readBedFile("/media/bs674/2t/openChromotainAndTFBS/all_reproducible_peaks_summits_merged.bed/all_reproducible_peaks_summits_merged.bed", tFBSBedRecords);
			
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_HiC_B73_leaf_merged_loops.txt", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/RNAPII.loop", hi_cRecords);
			readDnaLoopFile("/media/bs674/1_8t/AndCns/maizeHicAndMethylation/h3k4me3.loop", hi_cRecords);
			
//			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_husk.bed", h3k9acs);
//			readBedFile("/media/bs674/1_8t/AndCns/overlapWithH3k9/GSE94251_H3K9ac_ist.bed", h3k9acs);
			
			readBamFile(samFile, cnsRecords);
			System.out.println("bed files reading done");
			
			for( String chr : openChromatinRecords.keySet() ) {
				Collections.sort(tFBSBedRecords.get(chr));
				Collections.sort(openChromatinRecords.get(chr));
				Collections.sort(cnsRecords.get(chr));
				Collections.sort(h3k9acs.get(chr));
				Collections.sort(hi_cRecords.get(chr));
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					cnsRecord.setTag( (short)1 );
				}
//				for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
//					openChromatin.setTag( (short)2 );
//				}
//				for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
//					tFBSBedRecord.setTag( (short) 4 );
//				}
//				for( BedRecord hi_cRecord : hi_cRecords.get(chr) ) {
//					hi_cRecord.setTag( (short) 8 );
//				}
//				for( BedRecord h3k9ac : h3k9acs.get(chr) ) {
//					h3k9ac.setTag( (short) 16 );
//				}
			}
			
			for ( String chr : cnsRecords.keySet() ) {
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
						if ( openChromatin.overlapWith(cnsRecord)){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)2) );
							//openChromatin.setTag( (short) (openChromatin.getTag() | (short)1) );
							break;
						}else if ( cnsRecord.getEnd()<= openChromatin.getStart() ) {
							break;
						}
					}
					
					for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
						if ( tFBSBedRecord.overlapWith(cnsRecord)){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)4) );
							//tFBSBedRecord.setTag( (short) (tFBSBedRecord.getTag() | (short)1) );
							break;
						}else if ( cnsRecord.getEnd()<= tFBSBedRecord.getStart() ) {
							break;
						}
					}
					
					for( BedRecord h3k9ac : h3k9acs.get(chr) ) {
						if ( h3k9ac.overlapWith(cnsRecord)){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)8) );
							//h3k9ac.setTag( (short) (h3k9ac.getTag() | (short)1) );
							break;
						}else if ( cnsRecord.getEnd()<= h3k9ac.getStart() ) {
							break;
						}
					}
					
					for( BedRecord hi_cRecord : hi_cRecords.get(chr) ) {
						if ( hi_cRecord.overlapWith(cnsRecord)){
							cnsRecord.setTag( (short) (cnsRecord.getTag() | (short)16) );
							//hi_cRecord.setTag( (short) (hi_cRecord.getTag() | (short)1) );
							break;
						}else if ( cnsRecord.getEnd()<= hi_cRecord.getStart() ) {
							break;
						}
					}
				}
			}
			/*
			for ( String chr : openChromatinRecords.keySet() ) {
				for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
					for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
						if ( tFBSBedRecord.overlapWith(openChromatin)){
							openChromatin.setTag( (short) (openChromatin.getTag() | (short)4) );
							tFBSBedRecord.setTag( (short) (tFBSBedRecord.getTag() | (short)2) );
							break;
						}else if ( openChromatin.getEnd()<= tFBSBedRecord.getStart() ) {
							break;
						}
					}
				}
			}*/
			for ( String chr : cnsRecords.keySet() ) {
				HashMap<Integer, HashSet<String>> sequences = new HashMap<Integer, HashSet<String>>();
				HashMap<Integer, HashMap<String, String>> sequencesfull = new HashMap<Integer, HashMap<String, String>>();
				
				for( BedRecord cnsRecord : cnsRecords.get(chr) ) {
					if( ! sequences.containsKey((int) cnsRecord.getTag()) ) {
						sequences.put((int) cnsRecord.getTag(), new HashSet<String>());
						sequencesfull.put((int) cnsRecord.getTag(), new HashMap<String, String>());
					}
					if( (cnsRecord.getEnd() -cnsRecord.getStart()+1)>=seq_length  ) {
						int start = cnsRecord.getStart();
						int end = start + seq_length -1;
						sequences.get((int) cnsRecord.getTag()).add(maizeGenome.get(chr).substring(start-1, end));
						sequencesfull.get((int) cnsRecord.getTag()).put(maizeGenome2.get(chr).substring(start-1, cnsRecord.getEnd()), chr+":"+start+":"+cnsRecord.getEnd()+"\t"+cnsRecord.getQuerySequence());
					}
				}
			/*
				for( BedRecord openChromatin : openChromatinRecords.get(chr) ) {
					if( ! sequences.containsKey((int) openChromatin.getTag()) ) {
						sequences.put((int) openChromatin.getTag(), new HashSet<String>());
						sequencesfull.put((int) openChromatin.getTag(), new HashMap<String, String>());
					}
					if( (openChromatin.getTag() & 1) == 0  ) { // this is not a CNS
						if( (openChromatin.getEnd() -openChromatin.getStart()+1)>=seq_length  ) {
							int start = openChromatin.getStart();
							int end = start + seq_length -1;
							sequences.get((int) openChromatin.getTag()).add(maizeGenome.get(chr).substring(start-1, end));
							sequencesfull.get((int) openChromatin.getTag()).put(maizeGenome2.get(chr).substring(start-1, openChromatin.getEnd()), chr+":"+start+":"+openChromatin.getEnd());
						}
					}
				}
				for( BedRecord tFBSBedRecord : tFBSBedRecords.get(chr) ) {
					if( ! sequences.containsKey((int) tFBSBedRecord.getTag()) ) {
						sequences.put((int) tFBSBedRecord.getTag(), new HashSet<String>());
						sequencesfull.put((int) tFBSBedRecord.getTag(), new HashMap<String, String>());
					}
					if( (tFBSBedRecord.getTag() & 1) == 0 && (tFBSBedRecord.getTag() & 2) == 0  ) { // this is not a CNS and not a openChromatin
						if( (tFBSBedRecord.getEnd() -tFBSBedRecord.getStart()+1)>=seq_length  ) {
							int start = tFBSBedRecord.getStart();
							int end = start + seq_length -1;
							sequences.get((int) tFBSBedRecord.getTag()).add(maizeGenome.get(chr).substring(start-1, end));
							sequencesfull.get((int) tFBSBedRecord.getTag()).put(maizeGenome2.get(chr).substring(start-1, tFBSBedRecord.getEnd()), chr+":"+start+":"+tFBSBedRecord.getEnd());
						}
					}
				}*/
				
				HashMap<Integer, PrintWriter> outPuts = new HashMap<Integer, PrintWriter>();
				HashMap<Integer, PrintWriter> outPutFulls = new HashMap<Integer, PrintWriter>();
				for ( int tag : sequences.keySet() ) {
					if ( ! outPuts.containsKey(tag)  ) {
						PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/kmerClusterCns/chr"+ chr +"_sequences_with_tag_"+tag);
						PrintWriter outPutFull = new PrintWriter("/media/bs674/1_8t/AndCns/kmerClusterCns/chr"+ chr +"_full_sequences_with_tag_"+tag);
						outPuts.put(tag, outPut);
						outPutFulls.put(tag, outPutFull);
					}
					for( String sequence : sequences.get(tag) ) {
						outPuts.get(tag).println(tag + "\t" + sequence);
					}
					for( String sequence : sequencesfull.get(tag).keySet() ) {
						outPutFulls.get(tag).println(sequencesfull.get(tag).get(sequence) + "\t" + sequence);
					}
				}
				for ( int tag : outPuts.keySet() ) {
					outPuts.get(tag).close();
					outPutFulls.get(tag).close();
				}
			}
			
			HashMap<String, PrintWriter> outPuts = new HashMap<String, PrintWriter>();
			File file = new File( "/media/bs674/1_8t/AndCns/kmerClusterCns/B73.structuralTEv2.fulllength.2018-09-19.gff3" );
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				if( tempString.charAt(0) != '#' ) {
					String[] arrOfStr = tempString.split("\\s+");
					if ( ! outPuts.containsKey(arrOfStr[0])  ) {
						PrintWriter outPut = new PrintWriter("/media/bs674/1_8t/AndCns/kmerClusterCns/chr"+ arrOfStr[0] +"_sequences_with_tag_TE");
						outPuts.put(arrOfStr[0], outPut);
					}
					int start = Integer.parseInt(arrOfStr[3]);
					int end = Integer.parseInt(arrOfStr[4]);
					/*if( (end - start + 1)<seq_length ) {
						int length = (end - start + 1);
						start = (seq_length + length)/2 + start;
						if( start < 1  ) {
							start = 1;
						}
						end = start + seq_length - 1;
						if( end > maizeGenome.get(arrOfStr[0]).length() ) {
							end = maizeGenome.get(arrOfStr[0]).length();
						}
					}else */if( (end - start + 1) > seq_length ) {
						end = start + seq_length - 1;
						outPuts.get(arrOfStr[0]).println("TE\t"+maizeGenome.get(arrOfStr[0]).substring(start-1, end));
					}
				}
			}
			reader.close();
			for ( String chr : outPuts.keySet() ) {
				outPuts.get(chr).close();
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
	            	BedRecord bedRecord = new BedRecord(start, end);
	            	bedRecord.setQuerySequence(samRecord.getReadString());
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
