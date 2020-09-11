package me.songbx.model;

public class Alignment {
	private String referenceAlignment;
	private String queryAlignment;
	public String getReferenceAlignment() {
		return referenceAlignment;
	}
	public void setReferenceAlignment(String referenceAlignment) {
		this.referenceAlignment = referenceAlignment;
	}
	public String getQueryAlignment() {
		return queryAlignment;
	}
	public void setQueryAlignment(String queryAlignment) {
		this.queryAlignment = queryAlignment;
	}
	public Alignment(String referenceAlignment, String queryAlignment) {
		super();
		this.referenceAlignment = referenceAlignment.toUpperCase();
		this.queryAlignment = queryAlignment.toUpperCase();
	}
	public int getIdenticalBp() {
		int identicals = 0;
		for ( int i=0; i<referenceAlignment.length(); ++i ) {
			if ( referenceAlignment.charAt(i) == queryAlignment.charAt(i) ) {
				++identicals;
			}
		}
		return identicals;
	}
	public int getLargestIdenticalBp() {
		int largestIdenticalBp = 0;
		int thisIdenticals = 0;
		for ( int i=0; i<referenceAlignment.length(); ++i ) {
			if ( referenceAlignment.charAt(i) == queryAlignment.charAt(i) ) {
				++thisIdenticals;
			}else {
				if( thisIdenticals > largestIdenticalBp) {
					largestIdenticalBp = thisIdenticals;
				}
				thisIdenticals=0;
			}
		}
		return largestIdenticalBp;
	}
}
