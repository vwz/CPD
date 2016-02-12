package common;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

public class MathUtil {
	
	public static int Sum(ArrayList<String> list)
	{
		int sum = 0;
		for(String item: list)
		{
			sum += Integer.parseInt(item);
		}
		return sum;
	}

	class ValueComparator implements Comparator<Integer> {

	    HashMap<Integer, Float> base;
	    public ValueComparator(HashMap<Integer, Float> base) {
	        this.base = base;
	    }

	    // Note: this comparator imposes orderings that are inconsistent with equals.    
	    public int compare(Integer a, Integer b) {
	        if (base.get(a) >= base.get(b)) {
	            return -1;
	        } else {
	            return 1;
	        } // returning 0 would merge keys
	    }
	}
	
	
	
	public List<Integer> mapSortDouble(HashMap<Integer, Float> map){

		List<Integer> terms = new ArrayList<Integer>();
		ValueComparator bvc =  new ValueComparator(map);
		TreeMap<Integer, Float> listData = new TreeMap<Integer, Float>(bvc);

		listData.putAll(map);
		
		Iterator iter = listData.entrySet().iterator(); 

		for(int j=0;j<listData.size();j++)
		{
			Entry entry = (Entry) iter.next(); 
			terms.add((Integer)entry.getKey());
		}
		return terms;
	}
	
	public String mapSortDouble(HashMap<Integer, Float> map, int count){

		String res = "";
		ValueComparator bvc =  new ValueComparator(map);
		TreeMap<Integer, Float> listData = new TreeMap<Integer, Float>(bvc);

		listData.putAll(map);
		
		Iterator iter = listData.entrySet().iterator(); 

		for(int j=0;j<count;j++)
		{
			Entry entry = (Entry) iter.next(); 
			res += ((Integer)entry.getKey()+":"+(float)entry.getValue()+"\t");
		}
		
		res = res.substring(0, res.length()-1);
		return res;
	}
}
