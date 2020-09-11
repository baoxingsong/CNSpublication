package me.songbx.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class V3SnpsIdToV4Position {
	HashMap<String, Integer> snpIdToPosition;
	public void readVcfFile( ArrayList<String> vcfFiles ){
		for( String vcfFile : vcfFiles) {
			try {
	        	File file = new File(vcfFile);
	    		BufferedReader reader = new BufferedReader(new FileReader(file));
	            String tempString = null;
				
				while ((tempString = reader.readLine()) != null) {
					if(tempString.startsWith("#")){
						
					}else{
						String[] arrOfStr = tempString.split("\\s+");
						int position = Integer.parseInt(arrOfStr[1]);
						String snpId = arrOfStr[2];
						//System.out.println("line 26 " + snpId + " " + position);
						snpIdToPosition.put(snpId, position);
					}
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        }
		}
	}
	public V3SnpsIdToV4Position() {
		snpIdToPosition = new HashMap<String, Integer>();
		ArrayList<String> vcfFiles = new ArrayList<String>();
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr1.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr2.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr3.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr4.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr5.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr6.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr7.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr8.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr9.vcf");
		vcfFiles.add("/media/bs674/pan_and_non_asse/282V4/uplifted_APGv4/agpv4_hmp3_sites/agpv4_chr10.vcf");
		this.readVcfFile(vcfFiles);
	}
	public HashMap<String, Integer> getSnpIdToPosition() {
		return snpIdToPosition;
	}
	public void setSnpIdToPosition(HashMap<String, Integer> snpIdToPosition) {
		this.snpIdToPosition = snpIdToPosition;
	}
}
