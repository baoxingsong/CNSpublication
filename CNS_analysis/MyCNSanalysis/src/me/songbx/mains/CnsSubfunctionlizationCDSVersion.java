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

import me.songbx.impl.ChromoSomeReadImpl;
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.AlignmentMatch;
import me.songbx.model.GeneSimple;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


// this one not finished yet, do not run it
public class CnsSubfunctionlizationCDSVersion {
	public static void main(String[] args) {		
		new CnsSubfunctionlizationCDSVersion( true, args[0]);
		new CnsSubfunctionlizationCDSVersion( false, args[0]);
	}
	public CnsSubfunctionlizationCDSVersion( boolean core, String syngeneFile ){
		
		HashMap<String, GeneSimple> maizeGenes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		HashMap<String, String> maizeGeneSequence = new ChromoSomeReadImpl("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/gene.fa").getChromoSomeHashMap();
		
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
		samfiles.add("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/sorghum.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/A1025_08May2019/A1025.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/1013Chrysopogonserrulatus/A1013.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/Miscanthus_sinensis/msinensis.bam");
		samfiles.add("/media/bs674/1_8t/AndCns/maskGenomeForGenomeAlignment/sugarcane_tareploid/sugarCane.bam");
		
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter("/media/bs674/1_8t/AndCns/CnsSubfunctionlizationCDSVersion/subfunctionlization_CDS_"+core+"_"+syngeneFile);
	        Pattern p0 = Pattern.compile("^(\\d+)[HS]");
	        Pattern p1 = Pattern.compile("(\\d+)[HS]$");
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
					ArrayList<HashMap<Integer /*maize position*/, HashSet<AlignmentMatch> /*sorghum position*/ >> alignedMaizeFragments1 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the first gene
					ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>> alignedMaizeFragments2 = new ArrayList<HashMap<Integer, HashSet<AlignmentMatch>>>(); // the second gene
					for (int i=0; i<2; ++i) {
						String maizeChr = maizeGenes.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).getChromeSomeName();
						for (String samfile : samfiles) {
							SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(samfile));
							
							SAMRecordIterator it = reader.iterator();
							while (it.hasNext()) {
					            final SAMRecord samRecord = it.next();
					            if (samRecord.getReadName().equalsIgnoreCase(sorghum_gene_to_maize_genes.get(sorghumGene).get(i))){
						            Matcher m0=p0.matcher(samRecord.getCigarString());
						            Matcher m1=p1.matcher(samRecord.getCigarString());
						            int maizeStart = 1;
						            int maizeEnd = maizeGeneSequence.get(sorghum_gene_to_maize_genes.get(sorghumGene).get(i)).length();
						            int sorghumStart = samRecord.getAlignmentStart();   // is the position of the target genome, not necessarily sorghum
						            
						            if ( samRecord.getReadNegativeStrandFlag() ) {
						            	if(m1.find()){
					            			maizeStart = maizeStart + Integer.parseInt(m1.group(1));
						            	}
						            	if(m0.find()){
						            		maizeEnd = maizeEnd - Integer.parseInt(m0.group(1));
						            	}
						            } else {
						            	if(m0.find()){
							            	maizeStart = maizeStart + Integer.parseInt(m0.group(1));
						            	}
						            	if(m1.find()){
						            		maizeEnd = maizeEnd - Integer.parseInt(m1.group(1));
						            	}
						            }
						            if( maizeChr.charAt(0)!='B' ) {
					            		HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment = new HashMap<Integer, HashSet<AlignmentMatch>>();
					            		int sorghumPosition = sorghumStart;
					            		if ( samRecord.getReadNegativeStrandFlag() && !samRecord.getCigarString().equalsIgnoreCase("*") ) {
				            				for ( int maizePosition = maizeEnd; maizePosition >= maizeStart; ) {
					            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
					            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
					            				while (matcher.find()) {
					            					int length = Integer.parseInt( matcher.group(1));
					            					char cigerLetter =  matcher.group(2).charAt(0);
					            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
					            						for( int il=0; il<length; ++il ) {
				            								if ( ! alignedMaizeFragment.containsKey(maizePosition) ) {
				            									alignedMaizeFragment.put(maizePosition, new HashSet< AlignmentMatch>());
				            								}
				            								alignedMaizeFragment.get(maizePosition).add(new AlignmentMatch(samfile, samRecord.getContig(), sorghumPosition));
					            							--maizePosition;
					            							++sorghumPosition;
					            						}
					            					} else if ( cigerLetter=='I' ) {
					            						maizePosition-=length;
					            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
					            						sorghumPosition+=length;
					            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
					            					
					            					}else {
					            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
					            					}
					            				}
					            			}
				            			} else if (!samRecord.getCigarString().equalsIgnoreCase("*")) {
				            				for ( int maizePosition = maizeStart; maizePosition<=maizeEnd; ) {
					            				Pattern pCigar = Pattern.compile("(\\d+)([MIDNSHP=XB])");
					            				Matcher matcher = pCigar.matcher(samRecord.getCigarString());
					            				while (matcher.find()) {
					            					int length = Integer.parseInt( matcher.group(1));
					            					char cigerLetter =  matcher.group(2).charAt(0);
					            					if( cigerLetter=='M' || cigerLetter=='=' || cigerLetter=='X' ) {
					            						for( int il=0; il<length; ++il ) {
				            								if ( ! alignedMaizeFragment.containsKey(maizePosition) ) {
				            									alignedMaizeFragment.put(maizePosition, new HashSet< AlignmentMatch>());
				            								}
				            								alignedMaizeFragment.get(maizePosition).add(new AlignmentMatch(samfile, samRecord.getContig(), sorghumPosition));
					            							++maizePosition;
					            							++sorghumPosition;
					            						}
					            					} else if ( cigerLetter=='I' ) {
					            						maizePosition+=length;
					            					} else if ( cigerLetter=='D' || cigerLetter=='N' ) {
					            						sorghumPosition+=length;
					            					} else if ( cigerLetter=='S' || cigerLetter=='H' ) {
					            					
					            					}else {
					            						System.err.println("here we could not deal with the cigar:" + cigerLetter +" well, please contact the developper for updating");
					            					}
					            				}
					            			}
				            			}
						            	if( 0 == i ) {
						            		alignedMaizeFragments1.add(alignedMaizeFragment);
						            	}else {
						            		alignedMaizeFragments2.add(alignedMaizeFragment);
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
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment : alignedMaizeFragments1 ) {
							for( int position1 : alignedMaizeFragment.keySet() ) {
								if( ! toremove.containsKey(position1)  ) {
									toremove.put(position1, new HashSet<String>());
								}
								for( AlignmentMatch a : alignedMaizeFragment.get(position1) ) {
									toremove.get(position1).add(a.getSpecies());
								}	
							}
						}
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment : alignedMaizeFragments1 ) {
							for( int position : toremove.keySet() ) {
								if ( toremove.get(position).size()<5 ) {
									if( alignedMaizeFragment.containsKey(position) ) {
										alignedMaizeFragment.remove(position);
									}
								}
							}
						}
						
						toremove = new HashMap<Integer, HashSet<String>>();
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment : alignedMaizeFragments2 ) {
							for( int position1 : alignedMaizeFragment.keySet() ) {
								if( ! toremove.containsKey(position1)  ) {
									toremove.put(position1, new HashSet<String>());
								}
								for( AlignmentMatch a : alignedMaizeFragment.get(position1) ) {
									toremove.get(position1).add(a.getSpecies());
								}
							}
						}
						for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment : alignedMaizeFragments2 ) {
							for( int position : toremove.keySet() ) {
								if ( toremove.get(position).size()<5 ) {
									if( alignedMaizeFragment.containsKey(position) ) {
										alignedMaizeFragment.remove(position);
									}
								}
							}
						}
					}
					
					HashSet<Integer> overLappedPositions = new HashSet<Integer>();
					HashSet<Integer> positions1 = new HashSet<Integer>();
					HashSet<Integer> positions2 = new HashSet<Integer>();
					for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment1 : alignedMaizeFragments1 ) {
						for( int position1 : alignedMaizeFragment1.keySet() ) {
							positions1.add(position1);
							outerloop:
							if( ! overLappedPositions.contains(position1) ) {
								for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment2 : alignedMaizeFragments2 ) {
									for( int position2 : alignedMaizeFragment2.keySet() ) {
										for( AlignmentMatch a : alignedMaizeFragment1.get(position1) ) {
											if( alignedMaizeFragment2.get(position2).contains(a) ) {
												overLappedPositions.add(position1);
												break outerloop;
											}
										}
									}
								}
							}
						}
					}
					
					for( HashMap<Integer, HashSet<AlignmentMatch>> alignedMaizeFragment : alignedMaizeFragments2 ) {
						for( int position : alignedMaizeFragment.keySet() ) {
							positions2.add(position);
						}
					}

					int this_all_position = (positions1.size() + positions2.size() - overLappedPositions.size() );
					outPut.println("CDS\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + this_all_position + "\t" + positions1.size() + "\t" + positions2.size());
//					System.out.println("CDS\t" + sorghumGene + "\t" + sorghum_gene_to_maize_genes.get(sorghumGene).get(0) + "\t" +sorghum_gene_to_maize_genes.get(sorghumGene).get(1) + "\t" + this_all_position + "\t" + positions1.size() + "\t" + positions2.size());
				}
			}
			outPut.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
	}
}

