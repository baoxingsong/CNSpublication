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
import me.songbx.impl.ReadCnsGenotype;
import me.songbx.impl.ReadExpressionRankFile;
import me.songbx.impl.ReadGffForSimpleGene;
import me.songbx.model.BedRecord;
import me.songbx.model.GeneSimple;
import me.songbx.model.Strand;
import me.songbx.service.MaizeV3GeneToV4AndGetPosition;

public class QuadraticAnalysis {
	public static void main( String argv[] ) {
		/*
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/groot_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/groot_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/groot_expression_derivation_from_B73_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/groot_expression_derivation_from_B73_V2", "/media/bs674/1_8t/AndCns/CNSBasedGwas/groot_expression_derivation_from_B73_quadraticAnalysis_V2_2k");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/gshoot_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/shoot_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/gshoot_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/gshoot_expression_derivation_from_B73_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/gshoot_expression_derivation_from_B73_V2", "/media/bs674/1_8t/AndCns/CNSBasedGwas/gshoot_expression_derivation_from_B73_quadraticAnalysis_V2_2k");
		*/
		/*
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_kern_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_kern_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_kern_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_Kern_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_kern_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_kern_expression_derivation_from_B73_quadraticAnalysis");
		*/
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_root_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_root_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_root_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_GRoot_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_root_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_GRoot_expression_derivation_from_B73_quadraticAnalysis");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_shoot_expression_rank","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_shoot_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_shoot_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_GShoot_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_shoot_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_shoot_expression_derivation_from_B73_quadraticAnalysis");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Base_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Base_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_l3base_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Base_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Base_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Base_expression_derivation_from_B73_quadraticAnalysis");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_expression_rank", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_expression_derivation_from_B73","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_L3Tip_expression_derivation_from_B73_quadraticAnalysis");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_expression_rank","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_expression","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAD_expression_derivation_from_B73_quadraticAnalysis");
		
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_expression_rank","/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_quadraticAnalysis");
		new QuadraticAnalysis("/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_expression_derivation_from_B73", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_expression", "/media/bs674/1_8t/AndCns/CNSBasedGwas/hai_LMAN_expression_derivation_from_B73_quadraticAnalysis");
		
	}
	public QuadraticAnalysis( String input, String input2, String output ) {
		
		HashSet<String> cns1 = new HashSet<String>();
		HashSet<String> cns2 = new HashSet<String>();
		HashSet<String> cns3 = new HashSet<String>();
		
		try {
			SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/CNSprofiles/nonIntronUtr_OpOrTFBS.sam"));
			SAMRecordIterator it = reader.iterator();
			while (it.hasNext()) {
				final SAMRecord samRecord = it.next();
	        	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
	        	cns1.add(seqName);
			}		
			reader.close();
			
			
			reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/CNSprofiles/loop.sam"));
			it = reader.iterator();
			while (it.hasNext()) {
				final SAMRecord samRecord = it.next();
	        	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
	        	cns2.add(seqName);
			}		
			reader.close();
			
			reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File("/media/bs674/1_8t/AndCns/CNSprofiles/unknownCNS.sam"));
			it = reader.iterator();
			while (it.hasNext()) {
				final SAMRecord samRecord = it.next();
	        	String seqName = samRecord.getContig() + ":" + samRecord.getStart() + "-" + samRecord.getEnd();
	        	cns3.add(seqName);
			}		
			reader.close();			
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		PrintWriter outPut = null;
		try {
			outPut = new PrintWriter(output);
			
			int distance = 2000; // for upstream 2kb
			
			ReadExpressionRankFile readExpressionRankFile = new ReadExpressionRankFile(input);
			HashMap<String/*gene*/, HashMap<String/*accession*/, Double/*expressionRank*/>> geneExpressionRank = readExpressionRankFile.getGeneExpressionRank();
			
			HashMap<String/*gene*/, HashMap<String/*accession*/, Double/*expressionRank*/>> geneExpression = new ReadExpressionRankFile(input2).getGeneExpressionRank();
			
			
			ReadCnsGenotype readCnsGenotype = new ReadCnsGenotype("/media/bs674/1_8t/AndCns/CNSBasedGwas/encodeCNSasGenotype/all.CNS.genotype.txt");
			HashMap<String/*chr*/, HashMap<String/* cns id*/, HashMap<String/*accession id*/, Integer /*CNS length in this accession*/>>> cnsGenotype = readCnsGenotype.getCnsGenotype();
			HashMap<String/*CNS id*/, BedRecord > cnsRecords = readCnsGenotype.getCnsRecords();
			HashMap<String/*CNS id*/, Integer > cnsLength = readCnsGenotype.getCnsLength();
			
			MaizeV3GeneToV4AndGetPosition maizeV3GeneToV4AndGetPosition = new MaizeV3GeneToV4AndGetPosition("/media/bs674/1_8t/AndCns/subfunctionlization/pancns/correlationCNSWithTissueSpecificExpression/v3_to_v4");
			HashMap<String, Integer> v3Counts = maizeV3GeneToV4AndGetPosition.getV3Counts();
			HashMap<String, String> v3ToV4 = maizeV3GeneToV4AndGetPosition.getV3ToV4();
			
			HashMap<String, GeneSimple> genes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
			
			for( String gene : geneExpressionRank.keySet() ) {
				if( v3Counts.containsKey(gene) && v3Counts.get(gene)==1 && genes.containsKey(v3ToV4.get(gene)) ) {
					String v4Gene = v3ToV4.get(gene);
					int start = genes.get(v4Gene).getStart();
					int end = genes.get(v4Gene).getEnd();
					String chr = genes.get(v4Gene).getChromeSomeName();
					Strand strand = genes.get(v4Gene).getStrand();
					ArrayList<String> cnsIds = new ArrayList<String>();
					ArrayList<String> cnsId1s = new ArrayList<String>();
					ArrayList<String> cnsId2s = new ArrayList<String>();
					ArrayList<String> cnsId3s = new ArrayList<String>();
					ArrayList<String> cnsId4s = new ArrayList<String>();
					if( cnsGenotype.containsKey(chr) ) {
						for( String cnsId : cnsGenotype.get(chr).keySet() ) {
							BedRecord cns = cnsRecords.get(cnsId);
							if( Strand.POSITIVE == strand ) {
								if( start>=cns.getEnd() && (start-distance)<cns.getEnd() ) {
									cnsIds.add(cnsId);
									if( cns1.contains(cnsId) ) {
										cnsId1s.add(cnsId);
									} else if( cns2.contains(cnsId) ) {
										cnsId2s.add(cnsId);
									} else if( cns3.contains(cnsId) ) {
										cnsId3s.add(cnsId);
									}else {
										cnsId4s.add(cnsId);
									}
								}
							}else {
								if( end<=cns.getStart() && (end+distance)>cns.getStart() ) {
									cnsIds.add(cnsId);
									if( cns1.contains(cnsId) ) {
										cnsId1s.add(cnsId);
									} else if( cns2.contains(cnsId) ) {
										cnsId2s.add(cnsId);
									} else if( cns3.contains(cnsId) ) {
										cnsId3s.add(cnsId);
									}else {
										cnsId4s.add(cnsId);
									}
								}
							}
						}
						int totalLength = 0;
						for( String cnsId : cnsIds ) {
							totalLength += cnsLength.get(cnsId);
						}
						
						int totalLength1 = 0;
						for( String cnsId : cnsId1s ) {
							totalLength1 += cnsLength.get(cnsId);
						}
						int totalLength2 = 0;
						for( String cnsId : cnsId2s ) {
							totalLength2 += cnsLength.get(cnsId);
						}
						int totalLength3 = 0;
						for( String cnsId : cnsId3s ) {
							totalLength3 += cnsLength.get(cnsId);
						}
						int totalLength4 = 0;
						for( String cnsId : cnsId4s ) {
							totalLength4 += cnsLength.get(cnsId);
						}
						
						for( String accession : geneExpressionRank.get(gene).keySet() ) {
							int thisTotalLength = 0;
							for( String cnsId : cnsIds ) {
								if( cnsGenotype.get(chr).get(cnsId).containsKey(accession) ) {
									thisTotalLength += cnsGenotype.get(chr).get(cnsId).get(accession);
								}else {
									System.err.println(accession);
								}
							}
							
							int thisTotalLength1 = 0;
							for( String cnsId : cnsId1s ) {
								if( cnsGenotype.get(chr).get(cnsId).containsKey(accession) ) {
									thisTotalLength1 += cnsGenotype.get(chr).get(cnsId).get(accession);
								}else {
									System.err.println(accession);
								}
							}
							
							int thisTotalLength2 = 0;
							for( String cnsId : cnsId2s ) {
								if( cnsGenotype.get(chr).get(cnsId).containsKey(accession) ) {
									thisTotalLength2 += cnsGenotype.get(chr).get(cnsId).get(accession);
								}else {
									System.err.println(accession);
								}
							}
							
							int thisTotalLength3 = 0;
							for( String cnsId : cnsId3s ) {
								if( cnsGenotype.get(chr).get(cnsId).containsKey(accession) ) {
									thisTotalLength3 += cnsGenotype.get(chr).get(cnsId).get(accession);
								}else {
									System.err.println(accession);
								}
							}
							
							int thisTotalLength4 = 0;
							for( String cnsId : cnsId4s ) {
								if( cnsGenotype.get(chr).get(cnsId).containsKey(accession) ) {
									thisTotalLength4 += cnsGenotype.get(chr).get(cnsId).get(accession);
								}else {
									System.err.println(accession);
								}
							}
							if( geneExpression.containsKey(gene) && geneExpression.get(gene).containsKey(accession) ) {
								
							
								outPut.println(gene+"\t"+accession+"\t"+geneExpressionRank.get(gene).get(accession)
									+"\t"+totalLength+"\t"+thisTotalLength+"\t"+totalLength1+"\t"+thisTotalLength1+
									"\t"+totalLength2+"\t"+thisTotalLength2+"\t"+totalLength3+"\t"+
									thisTotalLength3+"\t"+totalLength4+"\t"+thisTotalLength+"\t"+geneExpression.get(gene).get(accession));
							}else {
								System.err.println("gene:" + gene + " accession:" + accession);
							}
						}
					}
				}
			}
			outPut.close();
		} catch (IOException e) {
            e.printStackTrace();
        }
	}
}
