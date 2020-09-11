package me.songbx.mains;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.AlignmentMatch;
import me.songbx.model.GeneSimple;
import me.songbx.model.Strand;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CnsSubfunctionlizationPanCnsVersion {
	public static void main(String[] args) {
//		new CnsSubfunctionlizationV2(1500);
//		new CnsSubfunctionlizationV2(1600);
//		new CnsSubfunctionlizationV2(1700);
//		new CnsSubfunctionlizationV2(1800);
//		new CnsSubfunctionlizationV2(1900);
		new CnsSubfunctionlizationPanCnsVersion(2000, false, args[0]);
//		new CnsSubfunctionlizationV2(2100);
//		new CnsSubfunctionlizationV2(2200);
//		new CnsSubfunctionlizationV2(2300);
//		new CnsSubfunctionlizationV2(2400);
//		new CnsSubfunctionlizationV2(2500);
	}
	public CnsSubfunctionlizationPanCnsVersion( int maizeIntervalEnd, boolean core, String syngeneFile ){
		int maizeIntervalStart = 1;
		
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, String> maize = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/masked_B73_v4_k20_46.fa").getChromoSomeHashMap();
		
		// get the syntenic paralogous relationship between maize and sorghum begin
		HashMap<String, ArrayList<String>> sorghum_gene_to_maize_genes = new HashMap<String, ArrayList<String>>();
		try {
        	File file = new File(syngeneFile);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
			while ((tempString = reader.readLine()) != null) {
				if(tempString.startsWith("#")){
					
				}else{
					String[] arrOfStr = tempString.split("\\s+");
					if ( ! sorghum_gene_to_maize_genes.containsKey(arrOfStr[0]) ) {
						sorghum_gene_to_maize_genes.put(arrOfStr[0], new ArrayList<String>());
					}
					sorghum_gene_to_maize_genes.get(arrOfStr[0]).add(arrOfStr[1]);
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
		
		ArrayList<String> samfiles = new ArrayList<String> ();
		samfiles.add("/media/bs674/1_8t/AndCns/sorghum/and_cns_sorghum_maize_V2_smaller_kmer_frequency/5.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/A1025_08May2019/result/5.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/1013Chrysopogonserrulatus/result/5.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/Miscanthus_sinensis/result/5.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/sugarcane_tareploid/result/5.bam");
		
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/subfunctionlization/pancns/subfunctionlization_pan_pancns_"+maizeIntervalStart+"_"+maizeIntervalEnd+"_"+syngeneFile);
	        Pattern p = Pattern.compile("^(\\d+)H");
			for ( String sorghumGene : sorghum_gene_to_maize_genes.keySet() ) {
				if( sorghum_gene_to_maize_genes.get(sorghumGene).size() == 2 ) {
					{
						int i=0;
						if ( ! sorghum_gene_to_maize_genes.containsKey(sorghumGene) ) {
							System.out.println("line 63:" + sorghumGene);
							continue;
						}
						if ( ! maizeGenes.containsKey(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)) ) {
							System.out.println("line 67:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
						
						i=1;
						if ( ! maizeGenes.containsKey(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)) ) {
							System.out.println("line 82:" + sorghum_gene_to_maize_genes.get(sorghumGene).get(i));
							continue;
						}
					}
					
					// this also count number of Ms with cigar string , similar with samtools depth command
					int all_position=0;
					int all_position1=0;
					int all_position2=0;
					{ // upstream, upstream is defined from start codon
						
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments1 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the first gene
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments2 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the second gene
						for (int i=0; i<2; ++i) {
							int maizeStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() == Strand.NEGTIVE ) {
								maizeStrand = 0; // negative strand
							}
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeIntervalEnd;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-maizeIntervalStart;
							if ( 0 == maizeStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalStart;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalEnd;
							}
							if( maizeStart < 1 ) {
								maizeStart = 1;
							}
							for (String samfile : samfiles) {
								SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(samfile));
								System.out.println(maizeChr + " " + maizeStart + " " + maizeEnd);
								SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
								while (it.hasNext()) {
						            final SAMRecord samRecord = it.next();
						            Matcher m=p.matcher(samRecord.getCigarString());
						            if(m.find()){
							            if ( samRecord.getReadNegativeStrandFlag() ) {
						            		
							            	int queryStart = Integer.parseInt(m.group(1));
							            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
							            	int objectStart = samRecord.getAlignmentStart();
						            		
							            	if( maizeChr.charAt(0)!='B' ) {
							            		HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment = new HashMap<Integer, HashSet<AlignmentMatch>>();
						            			if ( samRecord.getReadNegativeStrandFlag() ) {
						            				int objectPosition = objectStart;
						            				for ( int queryPosition = queryEnd; queryPosition >= queryStart; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if(  maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
							            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
							            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
							            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
							            								}
							            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
							            							}
							            							--queryPosition;
							            							++objectPosition;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						queryPosition-=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						objectPosition+=length;
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			} else {
						            				int objectPosition = objectStart;
						            				for ( int queryPosition = queryStart; queryPosition<=queryEnd; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
							            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
							            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
							            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
							            								}
							            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
							            							}
							            							++queryPosition;
							            							++objectPosition;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						queryPosition+=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						objectPosition+=length;
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}
						            			
								            	if( 0 == i ) {
								            		alignedSorghumFragments1.add(alignedSorghumFragment);
								            	}else {
								            		alignedSorghumFragments2.add(alignedSorghumFragment);
								            	}
							            	}
							            }
									}
						        }
								try {
									reader.close();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
						
						if( core ) {
							HashMap<Integer, HashSet<String>> toremove = new HashMap<Integer, HashSet<String>>();
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}	
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
						}
						
						HashSet<Integer> overLappedPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment1 : alignedSorghumFragments1 ) {
							for( int position1 : alignedSorghumFragment1.keySet() ) {
								positions1.add(position1);
								outerloop:
								if( ! overLappedPositions.contains(position1) ) {
									for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment2 : alignedSorghumFragments2 ) {
										for( int position2 : alignedSorghumFragment2.keySet() ) {
											for( AlignmentMatch a : alignedSorghumFragment1.get(position1) ) {
												if( alignedSorghumFragment2.get(position2).contains(a) ) {
													overLappedPositions.add(position1);
													break outerloop;
												}
											}
										}
									}
								}
							}
						}
						
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.keySet() ) {
								positions2.add(position);
							}
						}

						all_position1 += positions1.size();
						all_position2 += positions2.size();
						
						int this_all_position = (positions1.size() + positions2.size() - overLappedPositions.size() );
						all_position += this_all_position;
						outPut.println("upstream\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + this_all_position + "\t" + positions1.size() + "\t" + positions2.size());
					}
					
					{ // intron
						
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments1 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the first gene
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments2 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the second gene
						
						for (int i=0; i<2; ++i) {
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getNumberOfCds() > 1 ) {
								
								String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
								int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart();
								int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd();
								if( maizeStart < 1 ) {
									maizeStart = 1;
								}
								
								for (String samfile : samfiles) {
									SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(samfile));
									
									SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
									while (it.hasNext()) {
							            final SAMRecord samRecord = it.next();
							            
							            Matcher m=p.matcher(samRecord.getCigarString());
							            if(m.find()){
							            	int strand = 1;
								            if ( samRecord.getReadNegativeStrandFlag() ) {
								            	strand = 0;
								            }
								            if ( strand == 0 ) {
							            		
								            	int queryStart = Integer.parseInt(m.group(1));
								            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
								            	int objectStart = samRecord.getAlignmentStart();
								            	
								            	if( maizeChr.charAt(0)!='B' ) {
								            		HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment = new HashMap<Integer, HashSet<AlignmentMatch>>();
								            		if ( samRecord.getReadNegativeStrandFlag() ) {
							            				int objectPosition = objectStart;
							            				for ( int queryPosition = queryEnd; queryPosition >= queryStart; ) {
								            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
								            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
								            				while (matcher.find()) {
								            					int length = Integer.parseInt( matcher.group(1));
								            					char cigerLetter =  matcher.group(2).charAt(0);
								            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
								            						for( int il=0; il<length; ++il ) {
								            							if( maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
								            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
								            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
								            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
								            								}
								            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
								            							}
								            							--queryPosition;
								            							++objectPosition;
								            						}
								            					} else if ( cigerLetter=='I' ) {
								            						queryPosition-=length;
								            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
								            						objectPosition+=length;
								            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
								            					
								            					}else {
								            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
								            					}
								            				}
								            			}
							            			}else {
							            				int objectPosition = objectStart;
							            				for ( int queryPosition = queryStart; queryPosition<=queryEnd; ) {
								            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
								            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
								            				while (matcher.find()) {
								            					int length = Integer.parseInt( matcher.group(1));
								            					char cigerLetter =  matcher.group(2).charAt(0);
								            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
								            						for( int il=0; il<length; ++il ) {
								            							if( maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
								            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
								            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
								            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
								            								}
								            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
								            							}
								            							++queryPosition;
								            							++objectPosition;
								            						}
								            					} else if ( cigerLetter=='I' ) {
								            						queryPosition+=length;
								            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
								            						objectPosition+=length;
								            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
								            					
								            					}else {
								            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
								            					}
								            				}
								            			}
							            			}
									            	if( 0 == i ) {
									            		alignedSorghumFragments1.add(alignedSorghumFragment);
									            	}else {
									            		alignedSorghumFragments2.add(alignedSorghumFragment);
									            	}
								            	}
								            }
										}
							        }
									try {
										reader.close();
									} catch (IOException e) {
										e.printStackTrace();
									}
								}
							}
						}
						if( core ) {
							HashMap<Integer, HashSet<String>> toremove = new HashMap<Integer, HashSet<String>>();
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}	
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
						}
						
						HashSet<Integer> overLappedPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment1 : alignedSorghumFragments1 ) {
							for( int position1 : alignedSorghumFragment1.keySet() ) {
								positions1.add(position1);
								outerloop:
								if( ! overLappedPositions.contains(position1) ) {
									for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment2 : alignedSorghumFragments2 ) {
										for( int position2 : alignedSorghumFragment2.keySet() ) {
											for( AlignmentMatch a : alignedSorghumFragment1.get(position1) ) {
												if( alignedSorghumFragment2.get(position2).contains(a) ) {
													overLappedPositions.add(position1);
													break outerloop;
												}
											}
										}
									}
								}
							}
						}
						
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.keySet() ) {
								positions2.add(position);
							}
						}


						all_position1 += positions1.size();
						all_position2 += positions2.size();
						all_position += (positions1.size() + positions2.size() - overLappedPositions.size() );
						int this_all_position = (positions1.size() + positions2.size() - overLappedPositions.size() );
						outPut.println("intron\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + this_all_position + "\t" + positions1.size() + "\t" + positions2.size());
						
					}
					{ // down stream
						
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments1 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the first gene
						ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedSorghumFragments2 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the second gene
						for (int i=0; i<2; ++i) {
							int maizeStrand = 1;
							if( maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStrand() == Strand.NEGTIVE ) {
								maizeStrand = 0; // negative strand
							}
							
							String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
							int maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalStart;
							int maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getEnd()+maizeIntervalEnd;
							if( maizeStart < 1 ) {
								maizeStart = 1;
							}
							if ( 0 == maizeStrand ) {
								maizeStart = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart() - maizeIntervalEnd;
								maizeEnd = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getStart()-maizeIntervalStart;
							}
							for (String samfile : samfiles) {
								SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(samfile));
								
								SAMRecordIterator it = reader.queryOverlapping(maizeChr, maizeStart, maizeEnd);
								while (it.hasNext()) {
						            final SAMRecord samRecord = it.next();
						            
						            Matcher m=p.matcher(samRecord.getCigarString());
						            if(m.find() ){
						            	int strand = 1;
							            if ( samRecord.getReadNegativeStrandFlag() ) {
							            	strand = 0;
							            }
							            if ( strand == 0 ) {
						            		
							            	int queryStart = Integer.parseInt(m.group(1));
							            	int queryEnd = queryStart + samRecord.getReadString().length()-1;
							            	int objectStart = samRecord.getAlignmentStart();
						            		
							            	if( maizeChr.charAt(0)!='B' ) {
							            		HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment = new HashMap<Integer, HashSet<AlignmentMatch>>();
							            		if ( samRecord.getReadNegativeStrandFlag() ) {
						            				int objectPosition = objectStart;
						            				for ( int queryPosition = queryEnd; queryPosition >= queryStart; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
							            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
							            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
							            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
							            								}
							            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
							            							}
							            							--queryPosition;
							            							++objectPosition;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						queryPosition-=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						objectPosition+=length;
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}else {
						            				int objectPosition = objectStart;
						            				for ( int queryPosition = queryStart; queryPosition<=queryEnd; ) {
							            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
							            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
							            				while (matcher.find()) {
							            					int length = Integer.parseInt( matcher.group(1));
							            					char cigerLetter =  matcher.group(2).charAt(0);
							            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
							            						for( int il=0; il<length; ++il ) {
							            							if( maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'N' &&
							            									maize.get(samRecord.getContig()).charAt(objectPosition-1)!= 'n' ) {
							            								if ( ! alignedSorghumFragment.containsKey(objectPosition) ) {
							            									alignedSorghumFragment.put(objectPosition, new HashSet< AlignmentMatch>());
							            								}
							            								alignedSorghumFragment.get(objectPosition).add(new AlignmentMatch(samfile, samRecord.getReadName(), queryPosition));
							            							}
							            							++queryPosition;
							            							++objectPosition;
							            						}
							            					} else if ( cigerLetter=='I' ) {
							            						queryPosition+=length;
							            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
							            						objectPosition+=length;
							            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
							            					
							            					}else {
							            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
							            					}
							            				}
							            			}
						            			}
						            			
								            	if( 0 == i ) {
								            		alignedSorghumFragments1.add(alignedSorghumFragment);
								            	}else {
								            		alignedSorghumFragments2.add(alignedSorghumFragment);
								            	}
							            	}
							            }
									}
						        }
								try {
									reader.close();
								} catch (IOException e) {
									e.printStackTrace();
								}
							}
						}
						if( core ) {
							HashMap<Integer, HashSet<String>> toremove = new HashMap<Integer, HashSet<String>>();
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}	
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position1 : alignedSorghumFragment.keySet() ) {
									if( ! toremove.containsKey(position1)  ) {
										toremove.put(position1, new HashSet<String>());
									}
									for( AlignmentMatch a : alignedSorghumFragment.get(position1) ) {
										toremove.get(position1).add(a.getSpecies());
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments1 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
							for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
								for( int position : toremove.keySet() ) {
									if ( toremove.get(position).size()<5 ) {
										if( alignedSorghumFragment.containsKey(position) ) {
											alignedSorghumFragment.remove(position);
										}
									}
								}
							}
						}
						
						HashSet<Integer> overLappedPositions = new HashSet<Integer>();
						HashSet<Integer> positions1 = new HashSet<Integer>();
						HashSet<Integer> positions2 = new HashSet<Integer>();
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment1 : alignedSorghumFragments1 ) {
							for( int position1 : alignedSorghumFragment1.keySet() ) {
								positions1.add(position1);
								outerloop:
								if( ! overLappedPositions.contains(position1) ) {
									for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment2 : alignedSorghumFragments2 ) {
										for( int position2 : alignedSorghumFragment2.keySet() ) {
											for( AlignmentMatch a : alignedSorghumFragment1.get(position1) ) {
												if( alignedSorghumFragment2.get(position2).contains(a) ) {
													overLappedPositions.add(position1);
													break outerloop;
												}
											}
										}
									}
								}
							}
						}
						
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedSorghumFragment : alignedSorghumFragments2 ) {
							for( int position : alignedSorghumFragment.keySet() ) {
								positions2.add(position);
							}
						}

						all_position1 += positions1.size();
						all_position2 += positions2.size();
						all_position += (positions1.size() + positions2.size() - overLappedPositions.size() );
						int this_all_position = (positions1.size() + positions2.size() - overLappedPositions.size() );
						outPut.println("downstream\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + this_all_position + "\t" + positions1.size() + "\t" + positions2.size());
					
					}
					outPut.println("all\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + all_position + "\t" + all_position1 + "\t" + all_position2);
				}
			}
			outPut.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
	}
	public static ArrayList<String> cirgarToAlignment( String referenceSequence, int referenceBpPos, String querySequence, int queryBpPos, String cigar  ){
		String queryAlignment="";
		String referenceAlignment="";
		Pattern p = Pattern.compile("(\\d+)([MIDNSHP=XB])");
		Matcher matcher = p.matcher(cigar);
		while (matcher.find()) {
			int length = Integer.parseInt( matcher.group(1));
			char cigerLetter =  matcher.group(2).charAt(0);
			if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
				for( int il=0; il<length; ++il ) {
					queryAlignment += querySequence.substring(queryBpPos-1, queryBpPos);
					referenceAlignment += referenceSequence.substring(referenceBpPos-1, referenceBpPos);
					++queryBpPos;
					++referenceBpPos;
				}
			} else if ( cigerLetter=='I' ) {
				referenceAlignment += StringUtils.repeat("-", length);
				for( int il=0; il<length; ++il ) {
					referenceAlignment += StringUtils.repeat("-", length);
					queryAlignment += querySequence.substring(queryBpPos-1, queryBpPos);
					++queryBpPos;
				}
			} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
				referenceAlignment += referenceSequence.substring(referenceBpPos-1, referenceBpPos+length-1);
				queryAlignment += StringUtils.repeat("-", length);
				referenceBpPos+=length;
			} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
			}else {
				System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
			}
		}
		ArrayList<String> alignment = new ArrayList<String>();
		alignment.add(referenceAlignment);
		alignment.add(queryAlignment);
		return alignment;
	}
}


