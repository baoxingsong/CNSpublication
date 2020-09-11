package me.songbx.impl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import me.songbx.model.GeneSimple;
import me.songbx.model.Strand;

/**
 * The gene is defined as from start codon to stop codon
 * */


public class ReadGffForSimpleGene {
	public HashMap<String, GeneSimple> getGenes() {
		return genes;
	}

	public void setGenes(HashMap<String, GeneSimple> genes) {
		this.genes = genes;
	}

	private HashMap<String, GeneSimple> genes = new HashMap<String, GeneSimple>();
	
	public ReadGffForSimpleGene(String fileLocation) {
		HashMap<String, String> transcriptId_to_geneId_map = new HashMap<String, String>();
		try {
        	File file = new File(fileLocation);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            Pattern p = Pattern.compile("^(\\S*)\t(.*)\tmRNA\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(.*)ID=transcript:(\\S+?);Parent=(\\S+?);"); // for maize
           
			while ((tempString = reader.readLine()) != null) {
            	Matcher m=p.matcher(tempString);
				if(tempString.startsWith("#")){
					
				}else{
					Matcher mf = null;
					if(m.find()){
						mf = m;
					}
	                if(mf != null){
	                	String transcriptId = mf.group(9);
	                	transcriptId = transcriptId.replaceAll("=Chr\\d\\.", "");
	                	transcriptId = transcriptId.replaceAll("\\.mrna", "");
	                	transcriptId = transcriptId.replaceAll("=", "");
	                	
	                	String geneId = mf.group(10);
	                	geneId = geneId.replaceAll("gene:", "");
	                	transcriptId_to_geneId_map.put(transcriptId, geneId);
	                }
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
//		System.out.println("gff reading round one done");
        try {
        	File file = new File(fileLocation);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = null;
            Pattern p = Pattern.compile("^(\\S*)\t(.*)\tCDS\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(.*)Parent(.*?)[;,].*$");
			Pattern p1 = Pattern.compile("^(\\S*)\t(.*)\tCDS\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(.*)Parent(.*)$");//rice
			Pattern p2 = Pattern.compile("^(\\S*)\\s*(\\S*)\\s*CDS\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)PARENT=(\\S*)$");// for the annotation of Brassica rapa
			Pattern p3 = Pattern.compile("^(\\S*)\\s*(\\S*)\\s*CDS\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)\\s*(\\S*)gene_id \"ID=(\\S*)\"\\;.*$");// fot the annotation of fshgene http://mustang.biol.mcgill.ca:8885/download/S.irio/si_fgenesh.gtf
			Pattern p4 = Pattern.compile("^(\\S*)\t(.*)\tCDS\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(.*)transcript_id\\s*\"(\\w+?)\";");
			
			while ((tempString = reader.readLine()) != null) {
            	Matcher m=p.matcher(tempString);
				Matcher m_1=p1.matcher(tempString);
				Matcher m_2=p2.matcher(tempString);
				Matcher m_3=p3.matcher(tempString);
				Matcher m_4=p4.matcher(tempString);
				if(tempString.startsWith("#")){
					
				}else{
					Matcher mf = null;
					if(m.find()){
						mf = m;
					}else if(m_1.find()){
						mf = m_1;
					}else if(m_2.find()){
						mf = m_2;
					}else if(m_3.find()){
						mf = m_3;
					}else if(m_4.find()){
						mf = m_4;
					}
	                if(mf != null){
	                	int start = Integer.parseInt(mf.group(3));
	                	int end = Integer.parseInt(mf.group(4));
	                	
	                	if(start>end){
	                		int temp=start;
	                		start = end;
	                		end = temp;
	                	}
	                	
	                	String ChrId;
	                	if(mf == m_2 || mf==m_3){
	                		ChrId = mf.group(1);
	                	}else{
	                		ChrId = mf.group(1);
	                	}
	                	
	                	Strand strand;
	                	if("-".equals(mf.group(6))){
	                		strand = Strand.NEGTIVE;
	                	}else{
	                		strand = Strand.POSITIVE;
	                	}
	                	
	                	String transcriptId = mf.group(9);
	                	transcriptId = transcriptId.replaceAll("=Chr\\d\\.", "");
	                	transcriptId = transcriptId.replaceAll("\\.mrna", "");
	                	transcriptId = transcriptId.replaceAll("=", "");
	                	transcriptId = transcriptId.replaceAll("transcript:", "");
	                	if ( transcriptId_to_geneId_map.containsKey(transcriptId) ) {
	                		if( ! genes.containsKey(transcriptId_to_geneId_map.get(transcriptId)) ) {
	                			GeneSimple gene = new GeneSimple(strand, ChrId);
	                			genes.put(transcriptId_to_geneId_map.get(transcriptId), gene);
	                		}
	                		if( start < genes.get(transcriptId_to_geneId_map.get(transcriptId)).getStart() ) {
	                			genes.get(transcriptId_to_geneId_map.get(transcriptId)).setStart(start);
	                		}
	                		if( end > genes.get(transcriptId_to_geneId_map.get(transcriptId)).getEnd() ) {
	                			genes.get(transcriptId_to_geneId_map.get(transcriptId)).setEnd(end);
	                		}
	                		genes.get(transcriptId_to_geneId_map.get(transcriptId)).setNumberOfCds(genes.get(transcriptId_to_geneId_map.get(transcriptId)).getNumberOfCds()+1);
	                	}else {
	                		System.err.println("parsing gff problem. Could not find the parent gene ID for " + transcriptId);
	                	}
	                }
				}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
	}
	public static void main(String argv[]) {
		HashMap<String, GeneSimple> genes = new ReadGffForSimpleGene("/media/bs674/2t/genomeSequence/maize/Zea_mays.AGPv4.34.gff3").getGenes();
		System.out.println("dione");
		for( String geneId : genes.keySet() ) {
			System.out.println(genes.get(geneId).getChromeSomeName() + "\t" + geneId + "\t" + genes.get(geneId).getStrand() + "\t" + genes.get(geneId).getStart() + "\t" + genes.get(geneId).getEnd());
		}
	}
}
