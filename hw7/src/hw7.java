import java.util.ArrayList;
import java.util.Scanner;

public class hw7 {

    public static void main(String [] args){
        PathFinder pf = new PathFinder();
        //Scanner reader = new Scanner()
        //System.out.println(args[0]);
        pf.readInput(args[0]);
        System.out.println();

        System.out.println("M-Path Distance: " + pf.dist2Dest(0,3,2,1));
        System.out.println("M-Path Distance: " + pf.dist2Dest(0,3,2,2));


        ArrayList<Integer> path = pf.path2Dest(0,3,2,2);
        if(path == null) System.out.println("No path to destination");
        else{
            for(int i = 0; i < path.size(); i++) System.out.printf("%5s", path.get(i));
            System.out.println();
        }

        double[] distances = pf.dist2All(0, 2);
        for(int i = 0; i < distances.length; i++) System.out.println("distance to " + i + " is " + distances[i]);
        System.out.println();

        //Number of M-Paths
        System.out.println("Number of shortest paths: " + pf.noOfMPaths2Dest(0, 0, 2));


        System.out.println("Cost of mR-Tree: "  + pf.mRtreeCostFromSource(0,2));

        int[] parents = pf.mRtreeFromSource(0,2);
        for(int i = 0 ; i < parents.length; i++) System.out.println("parent of " + i + " is " + parents[i]);


        //System.out.println("M-Path Distance: " + pf.path2Dest(0,3,2,1));





    }


}
