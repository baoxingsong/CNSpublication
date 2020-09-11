package me.songbx.mains;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.CigarToAlignment;
import me.songbx.model.Alignment;
import me.songbx.model.BedRecord;
import me.songbx.model.NearByGene;
import me.songbx.service.IfIntron;
import me.songbx.service.IfNoCodingRna;

public class CNSProfile {
	
	class BpPosition{
		private NearByGene nearDistance;
		private int element;
		public int getElement() {
			return element;
		}
		public void setElement(int element) {
			this.element = element;
		}
		public NearByGene getNearDistance() {
			return nearDistance;
		}
		public void setNearDistance(NearByGene nearDistance) {
			this.nearDistance = nearDistance;
		}
		public BpPosition( NearByGene nearDistance, int element) {
			super();
			this.nearDistance = nearDistance;
			this.element=element;
		}
	}
	
	public CNSProfile( ArrayList<String> bamFiles, String output, String output2){
		
		HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
		
		HashSet<String> validatedChrs = new HashSet<String>();
		validatedChrs.add("1");
		validatedChrs.add("2");
		validatedChrs.add("3");
		validatedChrs.add("4");
		validatedChrs.add("5");
		validatedChrs.add("6");
		validatedChrs.add("7");
		validatedChrs.add("8");
		validatedChrs.add("9");
		validatedChrs.add("10");
		IfIntron ifIntron = new IfIntron("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3");
		System.out.println("IfIntron reading done");
		IfNoCodingRna ifNoCodingRna = new IfNoCodingRna("/media/bs674/1_8t/AndCns/noncodingRNA/maizeb73nonCodingRna.gtf");
		System.out.println("IfNoCodingRna reading done");
		HashMap<String, HashSet<BedRecord>> allRegions = new HashMap<String, HashSet<BedRecord>>(); // this one was designed to get all the unique CNS regions and is not used any more
		HashMap<String, HashMap<Integer, BpPosition>> bpPositions = new HashMap<String, HashMap<Integer, BpPosition>>();
		
		try {
			PrintWriter outPut = new PrintWriter(output);
			PrintWriter outPut2 = new PrintWriter(output2);
			for( String bamFile : bamFiles) {
				SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				SAMRecordIterator it = reader.iterator();
				while (it.hasNext()) {
		            final SAMRecord samRecord = it.next();
		            if( validatedChrs.contains(samRecord.getContig()) ) {
			            BedRecord bedRecord = new BedRecord(samRecord.getStart(), samRecord.getEnd());
			            if( !allRegions.containsKey(samRecord.getContig()) ) {
			            	allRegions.put(samRecord.getContig(), new HashSet<BedRecord>());
			            }
			            //if( allRegions.get(samRecord.getContig()).contains(bedRecord) ) {
			            	
			            //}else {
			            	int classiFication = 0;
			            	for( int position=samRecord.getStart(); position<=samRecord.getEnd(); ++position ) {
			            		int element = ifIntron.getElement ( samRecord.getContig(), position );
			            		if( element ==1 || element==2 ) { // 1 is intron, 2 is CDS
			            			classiFication = 1;
			            			break;
			            		}
			            		if( element ==3 ) {
			            			classiFication = element;
			            		}
			            	}
			            	allRegions.get(samRecord.getContig()).add(bedRecord);
			            	Alignment alignment = CigarToAlignment.getAlignment(maizeGenome, samRecord);
				            int length = alignment.getIdenticalBp(); // the length here is defined as the number of identites
				            int largestIdenticalBp = alignment.getLargestIdenticalBp();
				            int queryLength = samRecord.getAlignmentEnd() - samRecord.getAlignmentStart() + 1;
			            	if(  classiFication == 3 ) {
			            		//3 is UTR, negative is 5' UTR // it is correct, you should think about is carefully if you want to change it
			            		NearByGene nearDistance = ifIntron.getMinidistanceWithCds ( samRecord.getContig(), samRecord.getStart() , samRecord.getEnd());		            		
		            			outPut.println(queryLength + "\t" + length + "\t" + largestIdenticalBp + "\t" + classiFication + "\t" + nearDistance.getMiniDistance() + "\t" + nearDistance.getGene());
			            	}else if ( classiFication == 0 ){ // 0 isnonOpTfbsEqtlCns_profile intergenic a positive distance is downstream of gene
			            		NearByGene nearDistance = ifIntron.getMinidistanceWithTranscript ( samRecord.getContig(), samRecord.getStart() , samRecord.getEnd());
			            		if( nearDistance.getMiniDistance()>100000 ) {
		            				System.out.println(length + "\t" + classiFication + "\t" + bamFile + "\t" + nearDistance.getMiniDistance() + "\t" + nearDistance.getGene() + "\t" + samRecord.getSAMString() );
		            			}
			            		if (classiFication == 0 && ifNoCodingRna.isNocodingRna(samRecord.getContig(), samRecord.getStart() , samRecord.getEnd()) ) {
			            			outPut.println(queryLength + "\t" + length + "\t" + largestIdenticalBp + "\t4\t" + nearDistance.getMiniDistance() + "\t" + nearDistance.getGene()); // no-coding RNA
			            		}else{
			            			outPut.println(queryLength + "\t" + length + "\t" + largestIdenticalBp + "\t" + classiFication + "\t" + nearDistance.getMiniDistance() + "\t" + nearDistance.getGene());
			            		}
			            	}else {
			            		NearByGene nearDistance = ifIntron.getMinidistanceWithTranscript ( samRecord.getContig(), samRecord.getStart() , samRecord.getEnd());
			            		outPut.println(queryLength + "\t" + length + "\t" + largestIdenticalBp + "\t" + classiFication + "\t0\t"+nearDistance.getGene());
			            	}
			            	int position=samRecord.getStart();
			            	for ( int index=0; index<alignment.getReferenceAlignment().length(); ++index ) {
			            		if ( alignment.getReferenceAlignment().charAt(index) == alignment.getQueryAlignment().charAt(index) ) {
			            			if( ! bpPositions.containsKey(samRecord.getContig()) ) {
			            				bpPositions.put(samRecord.getContig(), new HashMap<Integer, BpPosition>());
			            			}
			            			NearByGene nearDistance = new NearByGene(0, "NA");
			            			classiFication = 0;
			            			if( bpPositions.get(samRecord.getContig()).containsKey(position) ) {
			            				classiFication = bpPositions.get(samRecord.getContig()).get(position).getElement();
			            				nearDistance = bpPositions.get(samRecord.getContig()).get(position).getNearDistance();
			            			}else {
			            				int element = ifIntron.getElement ( samRecord.getContig(), position );
	        							if( element ==1 || element==2 ) { // intron
	    			            			classiFication = 1;
	    			            		}
	    			            		if( element ==3 ) {
	    			            			classiFication = element;
		    			            		nearDistance = ifIntron.getMinidistanceWithCds ( samRecord.getContig(), position);
	    			            		}
	    			            		if (classiFication == 0 && ifNoCodingRna.isNocodingRna(samRecord.getContig(), samRecord.getStart() , samRecord.getEnd()) ) {
	    			            			classiFication=4; // no-coding RNA
        			            		}
	    			            		if( classiFication == 0 ) {
		    			            		nearDistance = ifIntron.getMinidistanceWithTranscript ( samRecord.getContig(), position);
	    			            		}
	    			            		BpPosition bpPosition = new BpPosition(nearDistance, classiFication);
	    			            		bpPositions.get(samRecord.getContig()).put(position, bpPosition);
			            			}
			            			outPut2.println(samRecord.getContig() +"\t" +  position + "\t" + classiFication + "\t" + nearDistance.getMiniDistance());
			            		}
			            		if( alignment.getReferenceAlignment().charAt(index) !=  '-' ) {
			            			++position;
			            		}
	            			}
			            //}
					}
				}
				reader.close();
			}
			outPut.close();
			outPut2.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	public static void main(String[] args) {
//		ArrayList<String> bamFiles = new ArrayList<String>();
//		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.sam");
//		ArrayList<String> bamFiles = new ArrayList<String>();
//		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlCns.sam");
//		ArrayList<String> bamFiles = new ArrayList<String>();
//		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/sorghum_nonOpTfbsEqtlCns.sam");
//		
		
		ArrayList<String> bamFiles = new ArrayList<String>();
		bamFiles.add("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5.bam");
		bamFiles.add("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5.bam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/CNS_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/CNS_bp_profile");
		
		bamFiles.clear();
		bamFiles.add("/media/bs674/1_8t/AndCns/Setaria_italica/result/5.bam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/santeria_CNS_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/santeria_CNS_bp_profile");
		
		bamFiles.clear();
		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.sam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_CNS_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_CNS_bp_profile");
		
		bamFiles.clear();
		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlCns.sam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlCns_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlCns_bp_profile");
		
		bamFiles.clear();
		bamFiles.add("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.sam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns_bp_profile");
		bamFiles.clear();
		bamFiles.add("/media/bs674/1_8t/AndCns/Setaria_italica/and_cns_setaria_maize_V3/5.bam");
		new CNSProfile(bamFiles, "/media/bs674/1_8t/AndCns/CNSprofiles/santeria_v3_CNS_profile", "/media/bs674/1_8t/AndCns/CNSprofiles/santeria_v3_CNS_bp_profile");
		
	}
}
