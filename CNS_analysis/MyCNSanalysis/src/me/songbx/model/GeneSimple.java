package me.songbx.model;

public class GeneSimple implements Comparable<GeneSimple> {
	private Strand strand;
	private String chromeSomeName;
	private int start;
	private int end;
	private int numberOfCds;
	private String name;
	
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public GeneSimple(Strand strand, String chromeSomeName, int start, int end) {
		super();
		this.strand = strand;
		this.chromeSomeName = chromeSomeName;
		this.start = start;
		this.end = end;
		this.numberOfCds=0;
	}
	
	public GeneSimple(Strand strand, String chromeSomeName, int start, int end, String name) {
		super();
		this.strand = strand;
		this.chromeSomeName = chromeSomeName;
		this.start = start;
		this.end = end;
		this.name=name;
	}
	
	public GeneSimple(Strand strand, String chromeSomeName) {
		super();
		this.strand = strand;
		this.chromeSomeName = chromeSomeName;
		this.start = Integer.MAX_VALUE;
		this.end = Integer.MIN_VALUE;
		this.numberOfCds=0;
	}
	
	public int getNumberOfCds() {
		return numberOfCds;
	}

	public void setNumberOfCds(int numberOfCds) {
		this.numberOfCds = numberOfCds;
	}

	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	public String getChromeSomeName() {
		return chromeSomeName;
	}
	public void setChromeSomeName(String chromeSomeName) {
		this.chromeSomeName = chromeSomeName;
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
	@Override
	public int compareTo(GeneSimple arg0) {
		return this.getStart() - arg0.getStart();
	}
}
