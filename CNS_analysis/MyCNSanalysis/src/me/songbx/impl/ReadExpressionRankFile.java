package me.songbx.impl;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class ReadExpressionRankFile {
	//HashMap<String /*chr*/, ArrayList<BedRecord> > cndRanges = new HashMap<String /*chr*/, ArrayList<BedRecord> >();
	HashMap<String/*gene*/, HashMap<String/*accession*/, Double/*expressionRank*/>> geneExpressionRank = new HashMap<String, HashMap<String, Double>>();
	
	public ReadExpressionRankFile( String fileLocation ) {
		try {
        	File file = new File(fileLocation);
    		BufferedReader reader = new BufferedReader(new FileReader(file));
            String tempString = reader.readLine();
            String [] accessions = tempString.split("\t");
            while ((tempString = reader.readLine()) != null) {
            	String [] elements = tempString.split("\t");
//            	assert( accessions.length == (elements.length-1));
//            	System.out.println("accessions: " + accessions.length + " elements: " + elements.length);
            	geneExpressionRank.put(elements[0], new HashMap<String, Double>());
            	for( int ei = 0; ei<accessions.length; ++ei ) {
            		geneExpressionRank.get(elements[0]).put(accessions[ei], Double.parseDouble(elements[ei+1]));
            	}
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
	}
	
	
	public HashMap<String, HashMap<String, Double>> getGeneExpressionRank() {
		return geneExpressionRank;
	}

	public void setGeneExpressionRank(HashMap<String, HashMap<String, Double>> geneExpressionRank) {
		this.geneExpressionRank = geneExpressionRank;
	}

	public static void main( String argv[] ) {
		HashMap<String/*gene*/, HashMap<String/*accession*/, Double/*expressionRank*/>> geneExpressionRank = (new ReadExpressionRankFile("/home/bs674/cbsuem02/local/workdir/bs674/CNSgwas/groot_expression_rank")).getGeneExpressionRank();
		for( String gene : geneExpressionRank.keySet() ) {
			for( String accession : geneExpressionRank.get(gene).keySet() ) {
				System.out.println(gene + "\t" + accession + "\t" + geneExpressionRank.get(gene).get(accession) );
			}
		}
	}
}


