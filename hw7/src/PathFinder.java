import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class PathFinder {

    int numIntersections;
    int numRoads;
    ArrayList<Intersection> intersections;
    ArrayList<Road> roads;
    HashMap<Intersection, Intersection> edges;


    public void readInput(String filePath) {
        File input = new File(filePath);
        Scanner reader = null;


        try {
            reader = new Scanner(input);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        int line = 0;
        while (reader.hasNextLine()) {
            String thisLine = reader.nextLine();

            thisLine = thisLine.trim();
            if (thisLine.equals("")) {
                continue;
            }
            //System.out.println(thisLine + " line: " + line);
            if (line == 0) {

                numIntersections = Integer.parseInt(thisLine.split("\\s+")[0]);
                intersections = new ArrayList<Intersection>(numIntersections);
//                System.out.println("numIntersections: " + numIntersections + " intersections size: " + intersections.toArray().length);
                numRoads = Integer.parseInt(thisLine.split("\\s+")[1]);
                roads = new ArrayList<Road>(numRoads);
            } else if (line > 0 && line < numIntersections + 1) {
                int first = Integer.parseInt(thisLine.split("\\s+")[0]);
                int second = Integer.parseInt(thisLine.split("\\s+")[1]);
                int third = Integer.parseInt(thisLine.split("\\s+")[2]);
                boolean added = false;
                for (int i = 0; i < intersections.size(); i++) {
                    if (intersections.get(i).getIdentifier() > first) {
                        added = true;
                        intersections.add(i, new Intersection(first, second, third));
                        break;
                    }
                }
                if (added == false) {
                    intersections.add(new Intersection(first, second, third));
                }
            } else if (line > numIntersections) {
                roads.add(new Road(intersections.get(Integer.parseInt(thisLine.split("\\s+")[0])), intersections.get(Integer.parseInt(thisLine.split("\\s+")[1]))));
                intersections.get(Integer.parseInt(thisLine.split("\\s+")[0])).addNeighbor(intersections.get(Integer.parseInt(thisLine.split("\\s+")[1])));
                edges.put(intersections.get(Integer.parseInt(thisLine.split("\\s+")[0])), intersections.get(Integer.parseInt(thisLine.split("\\s+")[1])));
            }
            line++;
        }
    }

    public PathFinder() {
        intersections = new ArrayList<Intersection>();
        roads = new ArrayList<Road>();
        edges = new HashMap<Intersection, Intersection>();
    }

    public double dist2Dest(int source, int destination, int ary, int dijkstraPath) {
        HashMap<Integer, Double> dist = new HashMap<Integer, Double>(numIntersections);
        HashMap<Intersection, Intersection> pred = new HashMap<Intersection, Intersection>(numIntersections);
        HashSet<Integer> visited = new HashSet<Integer>();


        if (dijkstraPath == 1) {
            Iterator<Intersection> iter = intersections.iterator();

            while (iter.hasNext()) {
                Intersection temp = iter.next();
                dist.put(temp.getIdentifier(), Double.MAX_VALUE);
            }

            dist.put(intersections.get(source).getIdentifier(), 0.0);

            MinHeap mh = new MinHeap(ary);
            //add the source node to the minheap

            mh.add(source, intersections.get(source).getDistance(intersections.get(source)));
            while (mh.isEmpty() == false) {

                ArrayList<Integer> testing = mh.getHeap();

                KeyValue temp = mh.peek();
                Intersection a = intersections.get(temp.getKey());
                if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                    mh.extractMin();
                    visited.add(a.getIdentifier());
                    Iterator<Intersection> itera = a.getNeighbors().iterator();
                    while (itera.hasNext()) {
                        Intersection iteraVertex = itera.next();
                        if (visited.contains(iteraVertex.getIdentifier()) == false) {
                            double alternative = dist.get(a.getIdentifier()).doubleValue() + a.getDistance(iteraVertex);
                            if (alternative < dist.get(iteraVertex.getIdentifier())) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                            }
                        }
                    }
                }
            }

            if (dist.get(intersections.get(destination).getIdentifier()) == Double.MAX_VALUE) {
                return -1.0;
            }

            return dist.get(intersections.get(destination).getIdentifier());

        } else if (dijkstraPath == 2) {
            Iterator<Intersection> iter = intersections.iterator();

            while (iter.hasNext()) {
                Intersection temp = iter.next();
                dist.put(temp.getIdentifier(), Double.MAX_VALUE);
            }

            dist.put(intersections.get(source).getIdentifier(), intersections.get(source).getDistance(intersections.get(destination)));

            MinHeap mh = new MinHeap(ary);
            //add the source node to the minheap

            mh.add(source, intersections.get(source).getDistance(intersections.get(destination)));
            while (mh.isEmpty() == false) {

                KeyValue temp = mh.peek();
                Intersection a = intersections.get(temp.getKey());
                if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                    mh.extractMin();
                    visited.add(a.getIdentifier());
                    Iterator<Intersection> itera = a.getNeighbors().iterator();
                    while (itera.hasNext()) {
                        Intersection iteraVertex = itera.next();
                        if (visited.contains(iteraVertex.getIdentifier()) == false) {
                            double alternative = dist.get(a.getIdentifier()).doubleValue() + a.getDistance(iteraVertex) + iteraVertex.getDistance(intersections.get(destination)) - a.getDistance(intersections.get(destination));
                            if (alternative < dist.get(iteraVertex.getIdentifier())) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                            }
                        }
                    }
                }
            }

            //the destination is not reachable
            if (dist.get(intersections.get(destination).getIdentifier()) == Double.MAX_VALUE) {
                return -1.0;
            }
            return dist.get(intersections.get(destination).getIdentifier());
        }

        return -1;
    }

    public ArrayList<Integer> path2Dest(int source, int destination, int ary, int dijkstraPath) {
        HashMap<Integer, Double> dist = new HashMap<Integer, Double>(numIntersections);
        HashMap<Integer, Integer> pred = new HashMap<Integer, Integer>(numIntersections);
        HashSet<Integer> visited = new HashSet<Integer>();


        if (dijkstraPath == 1) {
            Iterator<Intersection> iter = intersections.iterator();

            while (iter.hasNext()) {
                Intersection temp = iter.next();
                dist.put(temp.getIdentifier(), Double.MAX_VALUE);
                pred.put(temp.getIdentifier(), -1);
            }

            dist.put(intersections.get(source).getIdentifier(), 0.0);

            MinHeap mh = new MinHeap(ary);
            //add the source node to the minheap

            mh.add(source, intersections.get(source).getDistance(intersections.get(source)));
            while (mh.isEmpty() == false) {

                ArrayList<Integer> testing = mh.getHeap();

                KeyValue temp = mh.peek();
                Intersection a = intersections.get(temp.getKey());
                if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                    mh.extractMin();
                    visited.add(a.getIdentifier());
                    Iterator<Intersection> itera = a.getNeighbors().iterator();
                    while (itera.hasNext()) {
                        Intersection iteraVertex = itera.next();
                        if (visited.contains(iteraVertex.getIdentifier()) == false) {
                            double alternative = dist.get(a.getIdentifier()) + getDist(iteraVertex, a);
                            if (alternative < dist.get(iteraVertex.getIdentifier())) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                                pred.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                            } else if (alternative == dist.get(iteraVertex.getIdentifier()) && a.getIdentifier() < iteraVertex.getIdentifier()) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                                pred.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                            }
                        }
                    }
                }
            }

            if (dist.get(destination) == Double.MAX_VALUE) {
                return null;
            }

            //System.out.println("pred.tostring" + pred.toString());
            ArrayList<Integer> backTrack = new ArrayList<Integer>();
            int cur = destination;
            backTrack.add(cur);
            for (int i = 0; i < pred.size(); i++) {
                if (cur == source) {
                    break;
                }
                backTrack.add(0, pred.get(cur));
                cur = pred.get(cur);
            }

            return backTrack;

        } else if (dijkstraPath == 2) {
            Iterator<Intersection> iter = intersections.iterator();

            while (iter.hasNext()) {
                Intersection temp = iter.next();
                dist.put(temp.getIdentifier(), Double.MAX_VALUE);
                pred.put(temp.getIdentifier(), -1);
            }

            dist.put(intersections.get(source).getIdentifier(), intersections.get(source).getDistance(intersections.get(destination)));

            MinHeap mh = new MinHeap(ary);
            //add the source node to the minheap

            mh.add(source, intersections.get(source).getDistance(intersections.get(destination)));
            while (mh.isEmpty() == false) {

                KeyValue temp = mh.peek();
                Intersection a = intersections.get(temp.getKey());
                if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                    mh.extractMin();
                    visited.add(a.getIdentifier());
                    Iterator<Intersection> itera = a.getNeighbors().iterator();
                    while (itera.hasNext()) {
                        Intersection iteraVertex = itera.next();
                        if (visited.contains(iteraVertex.getIdentifier()) == false) {
                            double alternative = dist.get(a.getIdentifier()).doubleValue() + a.getDistance(iteraVertex) + iteraVertex.getDistance(intersections.get(destination)) - a.getDistance(intersections.get(destination));
                            if (alternative < dist.get(iteraVertex.getIdentifier())) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                                pred.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                            } else if (alternative == dist.get(iteraVertex.getIdentifier()) && a.getIdentifier() < iteraVertex.getIdentifier()) {
                                dist.replace(iteraVertex.getIdentifier(), alternative);
                                mh.add(iteraVertex.getIdentifier(), alternative);
                                pred.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                            }
                        }
                    }
                }
            }

            if (dist.get(destination) == Double.MAX_VALUE) {
                return null;
            }

            ArrayList<Integer> backTrack = new ArrayList<Integer>();
            int cur = destination;
            backTrack.add(cur);
            for (int i = 0; i < pred.size(); i++) {
                if (cur == source) {
                    break;
                }
                backTrack.add(0, pred.get(cur));
                cur = pred.get(cur);
            }

            return backTrack;
        }

        return null;
    }

    public double[] dist2All(int source, int ary) {
        HashMap<Integer, Double> dist = new HashMap<Integer, Double>(numIntersections);
        HashMap<Intersection, Intersection> pred = new HashMap<Intersection, Intersection>(numIntersections);
        HashSet<Integer> visited = new HashSet<Integer>();


        Iterator<Intersection> iter = intersections.iterator();

        while (iter.hasNext()) {
            Intersection temp = iter.next();
            dist.put(temp.getIdentifier(), Double.MAX_VALUE);
        }

        dist.put(intersections.get(source).getIdentifier(), 0.0);

        MinHeap mh = new MinHeap(ary);
        //add the source node to the minheap

        mh.add(source, intersections.get(source).getDistance(intersections.get(source)));
        while (mh.isEmpty() == false) {

            ArrayList<Integer> testing = mh.getHeap();

            KeyValue temp = mh.peek();
            Intersection a = intersections.get(temp.getKey());
            if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                mh.extractMin();
                visited.add(a.getIdentifier());
                Iterator<Intersection> itera = a.getNeighbors().iterator();
                while (itera.hasNext()) {
                    Intersection iteraVertex = itera.next();
                    if (visited.contains(iteraVertex.getIdentifier()) == false) {
                        double alternative = dist.get(a.getIdentifier()).doubleValue() + a.getDistance(iteraVertex);
                        if (alternative < dist.get(iteraVertex.getIdentifier())) {
                            dist.replace(iteraVertex.getIdentifier(), alternative);
                            mh.add(iteraVertex.getIdentifier(), alternative);
                        }
                    }
                }
            }
        }

        double[] distances = new double[dist.size()];

        for(int i = 0; i < distances.length; i++){
            if(dist.get(i) == Double.MAX_VALUE){
                distances[i] = -1;
                continue;
            }
            distances[i] = dist.get(i);
        }

        return distances;
    }

    public int noOfMPaths2Dest(int source, int destination, int ary){
        HashMap<Integer, Double> dist = new HashMap<Integer, Double>(numIntersections);
        HashMap<Intersection, Intersection> pred = new HashMap<Intersection, Intersection>(numIntersections);
        HashSet<Integer> visited = new HashSet<Integer>();


        Iterator<Intersection> iter = intersections.iterator();

        while (iter.hasNext()) {
            Intersection temp = iter.next();
            dist.put(temp.getIdentifier(), Double.MAX_VALUE);
        }

        dist.put(intersections.get(source).getIdentifier(), 0.0);

        MinHeap mh = new MinHeap(ary);
        //add the source node to the minheap

        mh.add(source, intersections.get(source).getDistance(intersections.get(source)));
        while (mh.isEmpty() == false) {

            ArrayList<Integer> testing = mh.getHeap();

            KeyValue temp = mh.peek();
            Intersection a = intersections.get(temp.getKey());
            if (a != null && !(dist.get(a.getIdentifier()) == Double.MAX_VALUE)) {
                mh.extractMin();
                visited.add(a.getIdentifier());
                Iterator<Intersection> itera = a.getNeighbors().iterator();
                while (itera.hasNext()) {
                    Intersection iteraVertex = itera.next();
                    if (visited.contains(iteraVertex.getIdentifier()) == false) {
                        double alternative = dist.get(a.getIdentifier()).doubleValue() + a.getDistance(iteraVertex);
                        if (alternative < dist.get(iteraVertex.getIdentifier())) {
                            dist.replace(iteraVertex.getIdentifier(), alternative);
                            mh.add(iteraVertex.getIdentifier(), alternative);
                        }
                    }
                }
            }
        }

        double curMin = Double.MAX_VALUE;
        int counter = 0;

        for(int i = 0; i < dist.size(); i++){
            if(dist.get(i) < curMin){
                curMin = dist.get(i);
                counter = 1;
            }
            else if(curMin == dist.get(i)){
                counter++;
            }

        }

        return counter;
    }

    public double getDist(Intersection a, Intersection b) {
        double x1 = a.getX();
        double y1 = a.getY();
        double x2 = b.getX();
        double y2 = b.getY();

        double xMinus = Math.abs(x1 - x2);
        double yMinus = Math.abs(y1 - y2);

        return Math.sqrt(xMinus * xMinus + yMinus * yMinus);
    }

    public int[] mRtreeFromSource(int source, int ary){

        //Going to keep track of the costs of the node with it's current attachment to the tree
        HashMap<Integer, Double> costs = new HashMap<Integer, Double>();

        //Going to keep track of the nodes that nodes are connected to the graph
        HashMap<Integer, Integer> parents = new HashMap<Integer, Integer>();
        HashSet<Integer> visited = new HashSet<Integer>();
        //
        Iterator<Intersection> iter = intersections.iterator();

        while (iter.hasNext()) {
            Intersection temp = iter.next();
            costs.put(temp.getIdentifier(), Double.MAX_VALUE);
            parents.put(temp.getIdentifier(), -1);
        }

        costs.put(intersections.get(source).getIdentifier(), 0.0);
        parents.replace(source, source);
        MinHeap mh = new MinHeap(ary);
        mh.add(source, costs.get(source));

        while(mh.isEmpty() == false){
            KeyValue temp = mh.extractMin();
            Intersection a = intersections.get(temp.getKey());

            if(a != null && !(costs.get(a.getIdentifier()) == Double.MAX_VALUE)){
                visited.add(a.getIdentifier());
                Iterator<Intersection> itera = a.getNeighbors().iterator();
                while(itera.hasNext()) {
                    Intersection iteraVertex = itera.next();
                    if (visited.contains(iteraVertex.getIdentifier()) == false && a.getDistance(iteraVertex) < costs.get(iteraVertex.getIdentifier())){
                        costs.replace(iteraVertex.getIdentifier(), a.getDistance(iteraVertex));
                        parents.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                        mh.add(iteraVertex.getIdentifier(), costs.get(iteraVertex.getIdentifier()));
                    }
                }
            }

        }

        int[] finalParents = new int[parents.size()];

        for (int i = 0; i < finalParents.length; i++){
            finalParents[i] = (int) parents.get(i);
        }

        return finalParents;
    }

    public double mRtreeCostFromSource(int source, int ary){
        //Going to keep track of the costs of the node with it's current attachment to the tree
        HashMap<Integer, Double> costs = new HashMap<Integer, Double>();

        //Going to keep track of the nodes that nodes are connected to the graph
        HashMap<Integer, Integer> parents = new HashMap<Integer, Integer>();
        HashSet<Integer> visited = new HashSet<Integer>();
        //
        Iterator<Intersection> iter = intersections.iterator();

        while (iter.hasNext()) {
            Intersection temp = iter.next();
            costs.put(temp.getIdentifier(), Double.MAX_VALUE);
            parents.put(temp.getIdentifier(), -1);
        }

        costs.put(intersections.get(source).getIdentifier(), 0.0);
        parents.replace(source, source);
        MinHeap mh = new MinHeap(ary);
        mh.add(source, costs.get(source));

        while(mh.isEmpty() == false){
            KeyValue temp = mh.extractMin();
            Intersection a = intersections.get(temp.getKey());

            if(a != null && !(costs.get(a.getIdentifier()) == Double.MAX_VALUE)){
                visited.add(a.getIdentifier());
                Iterator<Intersection> itera = a.getNeighbors().iterator();
                while(itera.hasNext()) {
                    Intersection iteraVertex = itera.next();
                    if (visited.contains(iteraVertex.getIdentifier()) == false && a.getDistance(iteraVertex) < costs.get(iteraVertex.getIdentifier())){
                        costs.replace(iteraVertex.getIdentifier(), a.getDistance(iteraVertex));
                        parents.replace(iteraVertex.getIdentifier(), a.getIdentifier());
                        mh.add(iteraVertex.getIdentifier(), costs.get(iteraVertex.getIdentifier()));
                    }
                }
            }

        }

        double totalCost = 0;

        for (int i = 0; i < costs.size(); i++){
            if(costs.get(i) == Double.MAX_VALUE){
                continue;
            }
            totalCost += costs.get(i);
        }

        return totalCost;
    }

    /**
     * Class represents a key value pair for a MinHeap, might've been able to use a hashmap, but who knows
     */
    class KeyValue {
        int key;
        double value;

        public KeyValue(int key, double value) {
            this.key = key;
            this.value = value;
        }

        public int getKey() {
            return key;
        }

        public double getValue() {
            return value;
        }

        public void setKey(int key) {
            this.key = key;
        }

        public void setValue(double value) {
            this.value = value;
        }
    }

    /**
     * Class that effectively represents a MinHeap for the implementation in Dijkstra's and Prim's
     */
    class MinHeap {
        KeyValue[] backing;
        int ary;
        int numElements;


        /**
         * Default Constructor
         */
        public MinHeap() {
            ary = 2;
            backing = new KeyValue[10];
            numElements = -1;
        }

        /**
         * Constructs a k-ary MinHeap with k being the integer passed into the constructor
         * @param ary
         */
        public MinHeap(int ary) {
            this.ary = ary;
            backing = new KeyValue[10];
            numElements = -1;
        }

        public boolean isEmpty() {
            if (numElements == -1) {
                return true;
            }
            return false;
        }

        public KeyValue peek() {
            if (isEmpty()) {
                return null;
            }
            return backing[0];
        }

        /**
         * Adds the key value pair to this MinHeap
         * @param key
         * @param value
         */
        public void add(int key, double value) {
            numElements++;
            expand();
            backing[numElements] = new KeyValue(key, value);

            heapUp(numElements);
        }

        public void heapUp(int index) {
            int smallest = (index - 1) / ary;

            while (smallest >= 0) {
                if (backing[index].getValue() < backing[smallest].getValue()) {
                    swap(index, smallest);
                    index = smallest;
                    smallest = (index - 1) / ary;
                } else if (backing[index].getValue() == backing[smallest].getValue() && backing[index].getKey() < backing[smallest].getKey()) {
                    swap(index, smallest);
                    index = smallest;
                    smallest = (index - 1) / ary;
                } else
                    break;
            }


        }

        public void heapDown(int index) {

            //keeps track of the indexes of the children of this index
            int[] children = new int[ary + 1];

            while (true) {
                for (int i = 1; i <= ary; i++) {
                    if (ary * index + i < numElements) {
                        children[i] = ((ary * index) + i);
                    } else {
                        children[i] = -1;
                    }

                }

                double minChild = Double.MAX_VALUE;
                int minChildIndex = -1;

                for (int i = 1; i <= ary; i++) {
                    if (children[i] != -1 && backing[children[i]].getValue() < minChild) {
                        minChild = backing[children[i]].getValue();
                        minChildIndex = children[i];
                    } else if (children[i] != -1 && backing[children[i]].getValue() == minChild) {
                        if (backing[children[i]].getKey() < backing[minChildIndex].getKey()) {
                            minChild = backing[children[i]].getValue();
                            minChildIndex = children[i];
                        }
                    }
                }

                if (minChild == Double.MAX_VALUE) {
                    break;
                }

                if (backing[index].getValue() > backing[minChildIndex].getValue())
                    swap(index, minChildIndex);
                else if (backing[index].getValue() == backing[minChildIndex].getValue() && backing[index].getKey() > backing[minChildIndex].getKey()) {
                    swap(index, minChildIndex);
                }
                index = minChildIndex;
            }

        }

        public KeyValue extractMin() {
            if (numElements == -1) {
                System.out.println("heap underflow");
            }

            KeyValue min = backing[0];
            backing[0] = backing[numElements];
            numElements = numElements - 1;
            heapDown(0);
            contract();
            return min;
        }

        public void swap(int a, int b) {
            KeyValue temp = backing[a];
            backing[a] = backing[b];
            backing[b] = temp;
        }

        public void expand() {
            //the array is full
            if (numElements == backing.length - 1) {
                KeyValue[] temp = new KeyValue[backing.length * 2];
                for (int i = 0; i < backing.length; i++) {
                    temp[i] = backing[i];
                }
                backing = temp;
            }
        }

        public void contract(){
            //the array is 1/4 full
            if(numElements < backing.length / 4){
                KeyValue[] temp = new KeyValue[backing.length / 2];
                for(int i = 0; i < temp.length; i++){
                    temp[i] = backing[i];
                }
                backing = temp;
            }
        }


        /**
         *
         * @return This MinHeap in the form of an arraylist of integers, the integers being the keys
         */
        public ArrayList<Integer> getHeap() {
            ArrayList<Integer> keys = new ArrayList<Integer>();
            for (int i = 0; i < numElements; i++) {
                keys.add(backing[i].getKey());
            }
            return keys;
        }
    }

    /**
     * Class effectively represents a node in a graph.
     */
    class Intersection {
        //unique identifier of this intersection
        int identifier;
        //coordinates of this intersection
        int[] coords;

        //intersections that this neighbor can connect to
        ArrayList<Intersection> neighbors;

        //edges connected to this intersection
        ArrayList<Road> edges;

        public Intersection() {
            identifier = -1;
            coords = null;
            neighbors = new ArrayList<Intersection>();
            edges = new ArrayList<Road>();
        }

        public Intersection(int identifier, int x, int y) {
            this.identifier = identifier;
            coords = new int[2];
            coords[0] = x;
            coords[1] = y;
            neighbors = new ArrayList<Intersection>();
            edges = new ArrayList<Road>();
        }

        public double getDistance(Intersection dist) {
            Road temp = new Road(new Intersection(identifier, coords[0], coords[1]), dist);
            temp.calculateWeight();
            return temp.getWeight();
        }

        public int getIdentifier() {
            return identifier;
        }

        public void setIdentifier(int identifier) {
            this.identifier = identifier;
        }

        public ArrayList<Intersection> getNeighbors() {
            return neighbors;
        }

        public ArrayList<Road> getEdges() {
            return edges;
        }

        public void setEdges(ArrayList<Road> edges) {
            this.edges = edges;
        }

        /**
         * sets the neighbors of this node
         *
         * @param neighbors
         */
        public void setNeighbors(ArrayList<Intersection> neighbors) {
            this.neighbors = neighbors;
        }

        /**
         * Adds an intersection to the list of neighbors for this intersection
         * If the neighbor is
         *
         * @param neighbor
         */
        public void addNeighbor(Intersection neighbor) {
            //the intersection does not have the given intersection as a neighbor
            if (neighbors.contains(neighbor) == false) {
                neighbors.add(neighbor);
                edges.add(new Road(new Intersection(identifier, getX(), getY()), neighbor));
            }

            //if the neighbor does not have
            if (neighbor.getNeighbors().contains(new Intersection(identifier, getX(), getY())) == false) {
                neighbor.getNeighbors().add(new Intersection(identifier, getX(), getY()));
                neighbor.addEdge(new Road(new Intersection(identifier, getX(), getY()), neighbor));
            }
        }

        /**
         * adds the edge to the start vertex and end vertex that this edge is attached to
         *
         * @param edge given edge
         */
        public void addEdge(Road edge) {
            if (edges.contains(edge)) {
                return;
            }
            edges.add(edge);
        }

        /**
         * if this intersection has an edge that goes to
         *
         * @param start
         * @param end
         * @return
         */
        public Road getEdge(Intersection start, Intersection end) {
            for (int i = 0; i < edges.size(); i++) {
                if (edges.get(i).isRoad(start, end) != -1) {
                    return edges.get(i);
                }
            }
            return null;
        }

        public int getX() {
            if (coords == null) {
                System.out.println("coords was null in getX");
            }
            return coords[0];
        }

        public int getY() {
            if (coords == null) {
                System.out.println("coords was null in getY");
            }
            return coords[1];
        }

        public void setX(int x) {
            if (coords == null) {
                coords = new int[2];
            }
            coords[0] = x;
        }

        public void setY(int y) {
            if (coords == null) {
                coords = new int[2];
            }
            coords[1] = y;
        }

        public String toString() {
            return "identifier: " + getIdentifier() + " x: " + getX() + " y: " + getY();
        }

    }

    /**
     * Class effectively represents an edge in a graph.
     */
    class Road {
        //weight of this road
        double weight;
        //beginning vertice of this road
        Intersection start;
        //end vertice of this road
        Intersection end;

        public double isRoad(Intersection start, Intersection end) {
            if (this.start == start && this.end == end) {
                return weight;
            } else if (this.end == start && this.start == end) {
                return weight;
            }
            return -1;
        }

        public Road() {
            weight = -1;
            start = null;
            end = null;
        }

        public Road(double weight, Intersection start, Intersection end) {
            this.weight = weight;
            this.start = start;
            this.end = end;
        }

        public Road(Intersection start, Intersection end) {
            calculateWeight();
            this.start = start;
            this.end = end;
        }

        /**
         * calculate the weight for this edge based on the start and end of this edge
         */
        public void calculateWeight() {
            if (start == null || end == null) {
                weight = -1;
                return;
            }

            weight = Math.sqrt((start.getX() - end.getX()) * (start.getX() - end.getX()) + (start.getY() - end.getY()) * (start.getY() - end.getY()));
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }

        public double getWeight() {
            return weight;
        }

        public void setStart(Intersection start) {
            this.start = start;
        }

        public Intersection getStart() {
            return start;
        }

        public void setEnd(Intersection end) {
            this.end = end;
        }

        public Intersection getEnd() {
            return end;
        }

        public String toString() {
            return "start: " + start.getIdentifier() + " end: " + end.getIdentifier();
        }
    }
}
