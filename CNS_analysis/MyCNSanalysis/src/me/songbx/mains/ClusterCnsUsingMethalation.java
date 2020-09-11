package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
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

public class ClusterCnsUsingMethalation {
	public static void main(String[] args) {
//		new ClusterCnsUsingMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.bam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns.methylationForEachRecord"
//				, "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_highMethylation.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtrCns_lowMethylation.sam");
		
//		new ClusterCnsUsingMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow.bam", "/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow.methylationForEachRecord"
//				, "/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow_highMethylation.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/stillUnknow_lowMethylation.sam");
		
		new ClusterCnsUsingMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.methylationForEachRecord");
		new ClusterCnsUsingMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/loop.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/loop.methylationForEachRecord");
		new ClusterCnsUsingMethalation("/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.sam", "/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.methylationForEachRecord");
		
	}
	
	public ClusterCnsUsingMethalation(String bamFile, String outputFile, String highMethylationSamOutput, String lowMethylationSamOutput ) {
		HashSet<String> mCG = new HashSet<String>();
		mCG.add("CGA");
		mCG.add("CGC");
		mCG.add("CGG");
		mCG.add("CGT");
		mCG.add("CGN");
		
		HashSet<String> CHG = new HashSet<String>();
		CHG.add("CAG");
		CHG.add("CCG");
		CHG.add("CTG");
		
		HashSet<String> CHH = new HashSet<String>();
		CHH.add("CAA");
		CHH.add("CCA");
		CHH.add("CTA");
		CHH.add("CAC");
		CHH.add("CCC");
		CHH.add("CTC");
		CHH.add("CAT");
		CHH.add("CCT");
		CHH.add("CTT");
		
		try {
			PrintWriter outPut = new PrintWriter(outputFile);
			PrintWriter outPut1 = new PrintWriter(highMethylationSamOutput);
			PrintWriter outPut2 = new PrintWriter(lowMethylationSamOutput);
			outPut.println("type\tseqName\tseqLength\tmCGs\tmCHGs\tmCHHs\tCGs\tCHGs\tCHHs");
			
			HashSet<String> unknowCns = new HashSet<String>();
			{
				SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/CNSprofiles/nonOpTfbsEqtlIntronUtrCns.bam"));
				SAMRecordIterator it = reader.iterator();
				while (it.hasNext()) {
		            final SAMRecord samRecord = it.next();
	            	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
	            	unknowCns.add(seqName);
				}
			}
			
			for( int chrInt=1; chrInt<=10; ++chrInt ){
				String chr = Integer.toString(chrInt);
				
				HashSet<Integer> mCGSites = new HashSet<Integer>();
				HashSet<Integer> mCHGSites = new HashSet<Integer>();
				HashSet<Integer> mCHHSites = new HashSet<Integer>();
				
				HashSet<Integer> mCGNoneSites = new HashSet<Integer>();
				HashSet<Integer> mCHGNoneSites = new HashSet<Integer>();
				HashSet<Integer> mCHHNoneSites = new HashSet<Integer>();
				
				File file = new File( "/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_allc_methylome_B73_leaf.tsv" );
				BufferedReader reader0 = new BufferedReader(new FileReader(file));
				String tempString = null;
				while ((tempString = reader0.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+");
					String methylationChr = arrOfStr[0];
					if( methylationChr.compareTo(chr) == 0 ) {
						String seq = arrOfStr[3];
						
						int position = Integer.parseInt( arrOfStr[1] );
						if( arrOfStr[6].charAt(0) == '1' ) {
							if( mCG.contains(seq) ) {
								mCGSites.add(position);
							}else if( CHG.contains(seq) ) {
								mCHGSites.add(position);
							}else if( CHH.contains(seq) ) {
								mCHHSites.add(position);
							}
						}else if( mCG.contains(seq) ) {
							mCGNoneSites.add(position);
						}else if( CHG.contains(seq) ) {
							mCHGNoneSites.add(position);
						}else if( CHH.contains(seq) ) {
							mCHHNoneSites.add(position);
						}
					}
				}
				reader0.close();
				System.out.println("mythylation " + chr + " reading done");
				
				
				SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				SAMRecordIterator it = reader.iterator();
				HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
				while (it.hasNext()) {
		            final SAMRecord samRecord = it.next();
		            if( samRecord.getContig().compareTo(chr) == 0 ) {
		            	Alignment alignment = CigarToAlignment.getAlignment(maizeGenome, samRecord);
		            	int seqLength = 0;
		            	int CGs = 0;
		            	int CHGs = 0;
		            	int CHHs = 0;
		            	int mCGs = 0;
		            	int mCHGs = 0;
		            	int mCHHs = 0;
		            	int chromoposition = samRecord.getAlignmentStart();
		            	boolean chromopositionChecked = false;
		            	for( int i=0; i<alignment.getReferenceAlignment().length(); ++i ) {
		            		if( alignment.getReferenceAlignment().charAt(i) == alignment.getQueryAlignment().charAt(i) ) {
		            			if( ! chromopositionChecked ) {
		            				if( mCGSites.contains(chromoposition) ) {
		            					++mCGs;
		            					++CGs;
			            			} else if ( mCGNoneSites.contains(chromoposition) ) {
			            				++CGs;
			            			} else if ( mCHGSites.contains(chromoposition) ) {
			            				++CHGs;
			            				++mCHGs;
			            			} else if ( mCHGNoneSites.contains(chromoposition) ) {
			            				++CHGs;
			            			} else if( mCHHSites.contains(chromoposition) ) {
			            				++CHHs;
			            				++mCHHs;
			            			} else if( mCHHNoneSites.contains(chromoposition) ) {
			            				++CHHs;
			            			}
		            				chromopositionChecked = true;
		            			}
		            			++seqLength;
		            		}
		            		if( alignment.getReferenceAlignment().charAt(i) != '-' ) {
		            			++chromoposition;
		            			chromopositionChecked = false;
		            		}
		            	}
		            	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
		            	if( unknowCns.contains(seqName) ) {
		            		outPut.print("unknown\t");
		            	} else {
		            		outPut.print("known\t");
		            	}
		            	outPut.println(seqName+"\t" + seqLength + "\t" + mCGs + "\t" + mCHGs + "\t" + mCHHs + "\t" + CGs + "\t" + CHGs + "\t" + CHHs);
		            	if( (((double) mCGs)/((double) CGs) > 0.5) || (((double) mCHGs)/((double) CHGs) > 0.5) ) {
		            		outPut1.print(samRecord.getSAMString());
		            	}else {
		            		outPut2.print(samRecord.getSAMString());
		            	}
		            }
				}	
			}
			outPut.close();
			outPut1.close();
			outPut2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public ClusterCnsUsingMethalation(String bamFile, String outputFile ) {
		HashSet<String> mCG = new HashSet<String>();
		mCG.add("CGA");
		mCG.add("CGC");
		mCG.add("CGG");
		mCG.add("CGT");
		mCG.add("CGN");
		
		HashSet<String> CHG = new HashSet<String>();
		CHG.add("CAG");
		CHG.add("CCG");
		CHG.add("CTG");
		
		HashSet<String> CHH = new HashSet<String>();
		CHH.add("CAA");
		CHH.add("CCA");
		CHH.add("CTA");
		CHH.add("CAC");
		CHH.add("CCC");
		CHH.add("CTC");
		CHH.add("CAT");
		CHH.add("CCT");
		CHH.add("CTT");
		
		try {
			PrintWriter outPut = new PrintWriter(outputFile);
			outPut.println("seqName\tseqLength\tmCGs\tmCHGs\tmCHHs\tCGs\tCHGs\tCHHs");
			
			for( int chrInt=1; chrInt<=10; ++chrInt ){
				String chr = Integer.toString(chrInt);
				
				HashSet<Integer> mCGSites = new HashSet<Integer>();
				HashSet<Integer> mCHGSites = new HashSet<Integer>();
				HashSet<Integer> mCHHSites = new HashSet<Integer>();
				
				HashSet<Integer> mCGNoneSites = new HashSet<Integer>();
				HashSet<Integer> mCHGNoneSites = new HashSet<Integer>();
				HashSet<Integer> mCHHNoneSites = new HashSet<Integer>();
				
				File file = new File( "/media/bs674/1_8t/AndCns/maizeHicAndMethylation/GSE120304_allc_methylome_B73_leaf.tsv" );
				BufferedReader reader0 = new BufferedReader(new FileReader(file));
				String tempString = null;
				while ((tempString = reader0.readLine()) != null) {
					String[] arrOfStr = tempString.split("\\s+");
					String methylationChr = arrOfStr[0];
					if( methylationChr.compareTo(chr) == 0 ) {
						String seq = arrOfStr[3];
						
						int position = Integer.parseInt( arrOfStr[1] );
						if( arrOfStr[6].charAt(0) == '1' ) {
							if( mCG.contains(seq) ) {
								mCGSites.add(position);
							}else if( CHG.contains(seq) ) {
								mCHGSites.add(position);
							}else if( CHH.contains(seq) ) {
								mCHHSites.add(position);
							}
						}else if( mCG.contains(seq) ) {
							mCGNoneSites.add(position);
						}else if( CHG.contains(seq) ) {
							mCHGNoneSites.add(position);
						}else if( CHH.contains(seq) ) {
							mCHHNoneSites.add(position);
						}
					}
				}
				reader0.close();
				System.out.println("mythylation " + chr + " reading done");
				
				
				SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
				SAMRecordIterator it = reader.iterator();
				HashMap<String, String> maizeGenome = new ChromoSomeReadImpl("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.dna.toplevel.fa").getChromoSomeHashMap();
				while (it.hasNext()) {
		            final SAMRecord samRecord = it.next();
		            if( samRecord.getContig().compareTo(chr) == 0 ) {
		            	Alignment alignment = CigarToAlignment.getAlignment(maizeGenome, samRecord);
		            	int seqLength = 0;
		            	int CGs = 0;
		            	int CHGs = 0;
		            	int CHHs = 0;
		            	int mCGs = 0;
		            	int mCHGs = 0;
		            	int mCHHs = 0;
		            	int chromoposition = samRecord.getAlignmentStart();
		            	boolean chromopositionChecked = false;
		            	for( int i=0; i<alignment.getReferenceAlignment().length(); ++i ) {
		            		if( alignment.getReferenceAlignment().charAt(i) == alignment.getQueryAlignment().charAt(i) ) {
		            			if( ! chromopositionChecked ) {
		            				if( mCGSites.contains(chromoposition) ) {
		            					++mCGs;
		            					++CGs;
			            			} else if ( mCGNoneSites.contains(chromoposition) ) {
			            				++CGs;
			            			} else if ( mCHGSites.contains(chromoposition) ) {
			            				++CHGs;
			            				++mCHGs;
			            			} else if ( mCHGNoneSites.contains(chromoposition) ) {
			            				++CHGs;
			            			} else if( mCHHSites.contains(chromoposition) ) {
			            				++CHHs;
			            				++mCHHs;
			            			} else if( mCHHNoneSites.contains(chromoposition) ) {
			            				++CHHs;
			            			}
		            				chromopositionChecked = true;
		            			}
		            			++seqLength;
		            		}
		            		if( alignment.getReferenceAlignment().charAt(i) != '-' ) {
		            			++chromoposition;
		            			chromopositionChecked = false;
		            		}
		            	}
		            	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
		            	outPut.println(seqName+"\t" + seqLength + "\t" + mCGs + "\t" + mCHGs + "\t" + mCHHs + "\t" + CGs + "\t" + CHGs + "\t" + CHHs);		            	
		            }
				}	
			}
			outPut.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
