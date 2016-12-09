import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Random;
import java.util.stream.Collectors;

/**
 * Generates a certain number of random vertices and forms a complete graph.
 * User has the ability to determine shortest path through Nearest Neighbor 
 * Algorithm or via a brute force method involving all permutations. For 
 * both shortest paths, the time of the algorithm along with the total
 * distance are shown.
 *   
 * @author VamsiG 
 * File: Graph.java
 */
public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices())
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  }



  // STUDENT CODE STARTS HERE

  /**
   * Forms n number of random vertices that are on a complete graph.
   * Edges are added on both sides to ensure that the graph is
   * undirected.
   * @param n the number of random vertices that are wanted 
   */
  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
    Random randomGenerator = new Random();

    //form all vertices
    for (int i = 0; i < n; i++ )
    {
       int x = randomGenerator.nextInt(100);
       int y = randomGenerator.nextInt(100);
       Vertex formedVertex = new Vertex(i, x, y);
       vertexNames.put(i, formedVertex);
    }

    //ensure completeness through adding edges between 
    //all points
    for (int i = 0; i < n; i++ )
    {
       Vertex currentVertex = vertexNames.get(i);

       for (int j = 0; j < n; j++ )
       {
         //if the objects aren't the same, make sure there's
         //an edge between them
         if (j != i)
         {
             Vertex target = vertexNames.get(j);
             Edge formedEdge = new Edge(currentVertex, target, 1);
             currentVertex.adjacentEdges.add(formedEdge);
         }
       }
    }
    
    computeAllEuclideanDistances(); // compute distances
  }

  /**
   * Runs the Nearest Neigbor algorithm on a graph. The initial 
   * point is chosen at random to allow for visualization of 
   * differences based upon starting vertex. 
   * @return List of edges that indicate the shortest route based
   * upon the nearest neighbor algorithm
   */
  public List<Edge> nearestNeighborTsp() {
    Random randomGenerator = new Random();
    //get random number between 0 and n-1 to 
    //get starting point for NN Algo
    int c = randomGenerator.nextInt(vertexNames.size() - 1);
    Vertex start = vertexNames.get(c);
    
    List<Edge> orderedEdges = new LinkedList<Edge>();

    //initialize all vertices
    for(Vertex v: vertexNames.values())
    {
      v.known = false;
    }

    //for all edges but the last edge (the edge from which the
    //last vertex returns to the intial edge), find the 
    //minimum edge on its adjacent edge list and then 
    //add that to the list indicating the route
    while(orderedEdges.size() < (vertexNames.values().size() - 1))
    {
      start.known = true;
      
      //determine the minimum distance and min index 
      //in the list of adjacent edges
      int minIndex = 0;
      double minDistance = 1000;
      for (int x = 0; x < start.adjacentEdges.size(); x++)
      {
        Edge check = start.adjacentEdges.get(x);
        if(check.target.known != true)
        {
          double checkDistance = check.distance;
        
          //if the new edge distance is smaller,
          //get its index
          if(checkDistance < minDistance)
          {
            minIndex = x;
            minDistance = checkDistance;
          }
        }
      }
      
      Edge minEdge = start.adjacentEdges.get(minIndex);
      //add min edge to list
      orderedEdges.add(minEdge);
      //change starting vertex for while loop
      start = minEdge.target;
    }

    //ensure that trip ends with return to vertex 0 by adding
    //that edge
    if(orderedEdges.size() == (vertexNames.values().size() - 1))
    {
      //get last vertex from existing list and 
      //add the edge that take it back to the initial 
      //point
      Vertex finalVertex = orderedEdges.get(orderedEdges.size() - 1).target;
      
      for(Edge e: finalVertex.adjacentEdges)
      {
        //make sure the edge is from final vertex to initial 
        //vertex
        if (e.target.name == orderedEdges.get(0).source.name)
        {
          orderedEdges.add(e);
        }
      }
    }
    
    return orderedEdges; // replace this line
  }

  /**
   * Runs a Brute Force algorithm on a graph to determine the 
   * shortest path possible. All paths are determined 
   * through permutations and then their distances are computed 
   * and compared to determine the minimum. 
   * @return List of edges that indicate the shortest route based
   * upon the brute force algorithm
   */
  public List<Edge> bruteForceTsp() 
  {
    
    LinkedList<int[]> permList = getPermutations(vertexNames.size());
    LinkedList<Double> distList = new LinkedList<Double>();
    List<Edge> shortestRoute = new LinkedList<Edge>();

    //for each permutation, send an array of permuted 
    //values to addEdges() and proceed to 
    //determine distances for each permutation
    for(int[] e: permList)
    {
      double sum = addEdges(e);
      distList.add(sum);
    }

    //determine the smallest permutation distance
    //and the index of its corresponding vertex
    //in the permutation list
    int minIndex = -1;
    double minValue = 1000000000;
    for (int x = 0; x < distList.size(); x++)
    {
        double check = distList.get(x);
      
        if(check < minValue)
        {
          minIndex = x;
          minValue = check;
        }
    }

    //form the list corressponding to the 
    //shortest route identified
    int[] minRoute = permList.get(minIndex);
    shortestRoute = formEdges(minRoute);

    return shortestRoute; 
  }

  /**
   * Takes a number and forms all permutations of 
   * the range to that number. This method utilizes 
   * recursion throuhg the detPerm method. 
   * @return Linked List of integer arrays containing all 
   * permutation of the path orders
   * @param n the range of ints that is required for the initial
   * array
   */
  public LinkedList<int[]> getPermutations(int n) 
  {
      LinkedList<int[]> permutedList = new LinkedList<int[]>();
      
      //form array corresponding to range from 0 to (n-1)
      int[] numArray = new int[n];
      for (int i = 0; i<n; i++)
      {
        numArray[i] = i;
      }

      //use detPerm() method for recursion to produce permutedList
      detPerm(numArray, 0, permutedList);
      return permutedList;
  }

  /**
   * Recursive method that determines all permutation. 
   * @param initArray initial array with values to be permuted
   * @param index current index
   * @param test LinkedList storing all the permuted arrays
   */
  private static void detPerm(int[] initArray, int index, LinkedList<int[]> test)
  {
      if(initArray.length - index == 1)
          //clone array to avoid the changes that will be made to it
          //java passes by values that represent references
          test.add(initArray.clone());
      else
          for(int i = index; i < initArray.length; i++){
              //swap index with changing i
              swap(initArray, index, i);
              //run algo again with index increased
              detPerm(initArray, index+1, test);
              //swap again after recursive call
              swap(initArray, index, i);
          }
  }

  /**
   * Swap elements at certain indeces within an array. 
   * @param arr initial array with values to be permuted
   * @param index1 first swap element's index
   * @param index2 second swap element's index
   */
  private static void swap(int[] arr, int index1, int index2)
  {
      int temp = arr[index1];
      arr[index1] = arr[index2];
      arr[index2] = temp;
  }

  /**
   * Form a list of edges based on an array of int values
   * @param n the array of ints representing vertices 
   * @return List of edges corresponding the the vertices 
   * reprented by the values in n
   */
  private List<Edge> formEdges (int[] n)
  {
    List<Edge> formedEdges = new LinkedList<Edge>();

    //deal with all edges but the last edge 
    //to return to initial point
    for(int x = 0; x < n.length - 1; x++)
    {
      int sourceVertValue = n[x];
      Vertex sourceVertex = vertexNames.get(sourceVertValue);

      int targetVertexValue = n[x + 1];

      //find the edge of source that corresponds with 
      //target value
      for(Edge m: sourceVertex.adjacentEdges)
      {
        if (m.target.name == targetVertexValue)
        {
          formedEdges.add(m);
          //break to avoid completion of for 
          //loop even after detection
          break;
        }
      }
    }

    //add edge from last point to initial point
    int sourceVertValue = n[n.length - 1];
    Vertex sourceVertex = vertexNames.get(sourceVertValue);

    int targetVertexValue = n[0];

    //find the edge of the last vertex that corresponds with 
    //the first vertex to complete loop
    for(Edge m: sourceVertex.adjacentEdges)
    {
      if (m.target.name == targetVertexValue)
      {
        formedEdges.add(m);
        //break to avoid completion of for 
        //loop even after detection
        break;
      }
    } 

    return formedEdges;
  }
    
  /**
   * Sum distances between edges in a list of integers
   * @param n the array of ints representing vertices 
   * @return the total distance between the vertices
   */
  private double addEdges (int[] n)
  {
    double sum = 0;

    //deal with all edges but the last edge 
    //to return to initial point
    for(int x = 0; x < n.length - 1; x++)
    {
      int sourceVertValue = n[x];
      Vertex sourceVertex = vertexNames.get(sourceVertValue);

      int targetVertexValue = n[x + 1];
      Vertex targetVertex = vertexNames.get(targetVertexValue);

      //add euclidean distance onto exiting sum
      sum = sum + computeEuclideanDistance(sourceVertex, targetVertex);
    }

    //add edge from last point to initial point
    int sourceVertValue = n[n.length - 1];
    Vertex sourceVertex = vertexNames.get(sourceVertValue);

    int targetVertexValue = n[0];
    Vertex targetVertex = vertexNames.get(targetVertexValue);

    //add euclidean distance onto exiting sum
    sum = sum + computeEuclideanDistance(sourceVertex, targetVertex);
    

    return sum;
  }

  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
