/**********************************************
  * Ryan Shea
  * MAE 423 - Heat Transfer
  * Final Project
**********************************************/

import java.util.ArrayList;

public class Flow
{
    private Node[][] field;
    private final double width;
    private final double height;
    private final double cellSize;
    private final int xNodes;
    private final int yNodes;
    private double flowTemperature = 300.0;
    private double surfaceTemperature = 400.0;
    private double reynoldsNumber = 200.0;
    private double kinematicViscosity = 0.00001569;
    private double alpha = 0.0002216;
    private double xLengthScale;
    public final double uInfinity;
    public final double uMax;
    public final double dt;
    
    private static final double EPSILON = 0.001;
    private static final int MAX_ITERATIONS = 1000;
    private static final double DT_SCALAR = 10;
    
    public enum BoundaryLocation {
        LEFT, RIGHT, BOTTOM, TOP
    }
    
    public enum NodeType {
        STREAM, INFLOW, OUTFLOW, FREELID, SOLIDWALL, SOLIDSURFACE, SOLIDINTERIOR//, POROUSWALL
    }

    private class Node
    {
        private NodeType type;
        private double psi; // streamfunction
        private double w; // vorticity
        private double u; // x velocity 
        private double v; // y velocity
        private double T; // temperature
        private double degreeOfSurfaceNormal = -1;
        
        public Node(int i, int j, NodeType type) {
            psi = uInfinity*j*cellSize;
            w = 0;
            u = 0;
            v = 0;
            T = flowTemperature;
            this.type = type;
        }
    }
    
    public Flow(double width, double height, double cellSize,
                double flowTemperature, double reynoldsNumber,
                double kinematicViscosity, double alpha, double xLengthScale)
    {
        this.width = width;
        this.height = height;
        this.cellSize = cellSize;
        this.xNodes = (int) Math.floor(width/cellSize) + 1;
        this.yNodes = (int) Math.floor(height/cellSize) + 1;
        this.xLengthScale = xLengthScale; // initially set as default
        this.field = new Node[xNodes][yNodes];
        
        this.flowTemperature = flowTemperature;
        this.reynoldsNumber = reynoldsNumber;
        this.kinematicViscosity = kinematicViscosity;
        this.alpha = alpha;
        this.uInfinity = reynoldsNumber*kinematicViscosity/xLengthScale;
        this.uMax = (10*kinematicViscosity/cellSize)*.75;
        this.dt = cellSize/(uMax*DT_SCALAR);
        
        createField(NodeType.INFLOW, NodeType.OUTFLOW, NodeType.FREELID, NodeType.FREELID);
    }
    
    public void createField(NodeType left, NodeType right,
                            NodeType bottom, NodeType top)
    {
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                if (i == 0) {
                    field[i][j] = new Node(i, j, left);
                } else if (i == xNodes-1) {
                    field[i][j] = new Node(i, j, right);
                } else if (j == 0) {
                    field[i][j] = new Node(i, j, bottom);
                } else if (j == yNodes-1) {
                    field[i][j] = new Node(i, j, top);
                } else {
                    field[i][j] = new Node(i, j, NodeType.STREAM);
                }
            }
        }
    }
    
    public void addSolidBox(double xCenter, double yCenter,
                            double halfWidth, double halfHeight,
                            double temperature)
    {
        xLengthScale = halfWidth*2;
        ArrayList<Node> interiorNodes = new ArrayList<Node>();
        double xMin = xCenter - halfWidth;
        double xMax = xCenter + halfWidth;
        double yMin = yCenter - halfHeight;
        double yMax = yCenter + halfHeight;
        for (int i = 1; i < xNodes-1; i++) {
            for (int j = 1; j < yNodes-1; j++) {
                if ((i*cellSize >= xMin && i*cellSize <= xMax) &&
                    (j*cellSize >= yMin && j*cellSize <= yMax)) {
                    field[i][j].type = NodeType.SOLIDSURFACE;
                    field[i][j].psi = uInfinity*yCenter;
                    field[i][j].T = temperature;
                    if ((i+1)*cellSize > xMax && i*cellSize <= xMax)
                        field[i][j].degreeOfSurfaceNormal = 0.0;
                    else if ((j+1)*cellSize > yMax && j*cellSize <= yMax)
                        field[i][j].degreeOfSurfaceNormal = 90.0;
                    else if ((i-1)*cellSize < xMin && i*cellSize >= xMin)
                        field[i][j].degreeOfSurfaceNormal = 180.0;
                    else if ((j-1)*cellSize < yMin && j*cellSize >= yMin)
                        field[i][j].degreeOfSurfaceNormal = 270.0;
                }
                if (field[i][j-1].type == NodeType.SOLIDSURFACE &&
                    field[i-2][j-1].type == NodeType.SOLIDSURFACE &&
                    field[i-1][j].type == NodeType.SOLIDSURFACE &&
                    field[i-1][j-2].type == NodeType.SOLIDSURFACE)
                {
                    interiorNodes.add(field[i-1][j-1]);
                }
            }
        }
        // clean up interior nodes
        for (Node node : interiorNodes) {
            node.type = NodeType.SOLIDINTERIOR;
            node.psi = 0;
        }
    }
    
    public void addSolidCylinder(double xCenter, double yCenter,
                                 double radius, double temperature)
    {
        xLengthScale = radius*2;
        double outerRadiusSquared = radius*radius;
        ArrayList<Node> interiorNodes = new ArrayList<Node>();
        for (int i = 1; i < xNodes-1; i++) {
            for (int j = 1; j < yNodes-1; j++) {
                double distanceToCenterSquared =
                    ((xCenter - i*cellSize)*(xCenter - i*cellSize) +
                     (yCenter - j*cellSize)*(yCenter - j*cellSize));
                if (distanceToCenterSquared <= outerRadiusSquared) {
                    field[i][j].type = NodeType.SOLIDSURFACE;
                    field[i][j].psi = uInfinity*yCenter;
                    field[i][j].T = temperature;
                    field[i][j].degreeOfSurfaceNormal =
                        degreeOnCircle(xCenter, yCenter, i*cellSize, j*cellSize);
                }
                if (field[i][j-1].type == NodeType.SOLIDSURFACE &&
                    field[i-2][j-1].type == NodeType.SOLIDSURFACE &&
                    field[i-1][j].type == NodeType.SOLIDSURFACE &&
                    field[i-1][j-2].type == NodeType.SOLIDSURFACE)
                {
                    interiorNodes.add(field[i-1][j-1]);
                }
            }
        }
        // clean up interior nodes
        for (Node node : interiorNodes) {
            node.type = NodeType.SOLIDINTERIOR;
            node.psi = 0;
        }
    }
    
    private void applyGaussSeidelIteration(double F)
    {
        double maxResidue;
        int iterations = 0;
        
        do {
            maxResidue = 0;
            for (int i = 0; i < xNodes; i++) {
                for (int j = 0; j < yNodes; j++) {
                    Node node = field[i][j];
                    switch (node.type) {
                        case STREAM:
                            double residueTerm =
                                (field[i+1][j].psi + field[i-1][j].psi +
                                 field[i][j+1].psi + field[i][j-1].psi +
                                 cellSize*cellSize*node.w -
                                 4*node.psi);
                            node.psi = node.psi + (F/4.0)*residueTerm;
                            if (Math.abs(residueTerm/node.psi) > maxResidue) {
                                maxResidue = Math.abs(residueTerm/node.psi);
                            }
                            break;
                        case OUTFLOW:
                            // update after convergence
                            break;
                        case INFLOW:
                        case FREELID:
                        case SOLIDWALL:
                        case SOLIDSURFACE:
                        case SOLIDINTERIOR:
                            // do nothing
                            break;
                        default:
                            System.out.println("Invalid node type in gauss seidel iteration");
                            System.exit(1);
                    }
                }
            }
            iterations++;
            if (iterations >= MAX_ITERATIONS) {
                System.out.println("Difficulty converging to gauss-seidel iterative solution.");
                break;
            }
        } while (maxResidue >= EPSILON);
        
        // update outflow
        for (int j = 0; j < yNodes; j++) {
            int i = xNodes-1;
            field[i][j].psi = 2*field[i-1][j].psi - field[i-2][j].psi;
        }
        
        updateVelocity();
    }
    
    private void updateVelocity()
    {
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                Node node = field[i][j];
                switch (node.type) {
                    case STREAM:
                    case OUTFLOW:
                    case FREELID:
                        if (j == 0) {
                            node.u = (field[i][j+1].psi - field[i][j].psi)/(cellSize);
                        } else if (j == yNodes-1) {
                            node.u = (field[i][j].psi - field[i][j-1].psi)/(cellSize);
                        } else {
                            node.u = (field[i][j+1].psi - field[i][j-1].psi)/(2*cellSize);
                        }
                        if (i == 0) {
                            node.v = -(field[i+1][j].psi - field[i][j].psi)/(cellSize);
                        } else if (i == xNodes-1) {
                            node.v = -(field[i][j].psi - field[i-1][j].psi)/(cellSize);
                        } else {
                            node.v = -(field[i+1][j].psi - field[i-1][j].psi)/(2*cellSize);
                        }
                        break;
                    case SOLIDWALL:
                    case SOLIDSURFACE:
                        node.v = 0;
                        node.u = 0;
                        break;
                    case INFLOW:
                    case SOLIDINTERIOR:
                        // do nothing
                        break;
                    default:
                        System.out.println("Invalid boundary condition in velocity calculation");
                        System.exit(1);
                }
            }
        }
        
        // update inflow
        for (int j = 0; j < yNodes; j++) {
            for (int i = 0; i < 3; i++) {
                field[i][j].u = uInfinity;
                field[i][j].v = 0;
            }
        }
    }
    
    private void advanceVorticity()
    {
        for (int i = 1; i < xNodes-1; i++) {
            for (int j = 1; j < yNodes-1; j++) {
                Node node = field[i][j];
                switch (node.type) {
                    case STREAM:
                        double delSquaredW = (field[i+1][j].w + field[i-1][j].w +
                                              field[i][j+1].w + field[i][j-1].w - 4*node.w)/(cellSize*cellSize);
                        double delUW;
                        double delVW;
                        if (node.u <= 0)
                            delUW = field[i+1][j].u*field[i+1][j].w - node.u*node.w;
                        else
                            delUW = node.u*node.w - field[i-1][j].u*field[i-1][j].w;
                        if (node.v <= 0)
                            delVW = field[i][j+1].v*field[i][j+1].w - node.v*node.w;
                        else
                            delVW = node.v*node.w - field[i][j-1].v*field[i][j-1].w;
                        node.w = node.w + dt*(-delUW/cellSize - delVW/cellSize +
                                              kinematicViscosity*delSquaredW);
                        break;
                    default:
                        break;
                }
            }
        }
        
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                Node node = field[i][j];
                switch (node.type) {
                    case SOLIDSURFACE:
                    case SOLIDWALL:
                        field[i][j].w = calculateStreamFunctionNormalDerivative(i, j);
                        break;
                    case OUTFLOW:
                        node.w = field[i-1][j].w;
                        break;
                    case INFLOW:
                    case FREELID:
                    case SOLIDINTERIOR:
                        node.w = 0;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    
    private static double degreeOnCircle(double centerX, double centerY,
                                        double p1X, double p1Y)
    {
        double dy = p1Y - centerY;
        double dx = p1X - centerX;
        double radians = Math.atan2(dy, dx);
        double degrees = Math.toDegrees(radians);
        if (degrees < 0) degrees += 360;
        return degrees;
    }
    
    private double calculateTemperatureNormalDerivative(int i, int j)
    {
        Node node = field[i][j];
        double degrees = node.degreeOfSurfaceNormal;
        if (i - 1 < 0 || j - 1 < 0 || i + 1 > xNodes-1 || j + 1 > xNodes-1)
            return 0;
        double d = getNormalValue(degrees, node.T,
                                  field[i-1][j].T, field[i+1][j].T,
                                  field[i][j-1].T, field[i][j+1].T);
        return 2*d/(cellSize*cellSize);
    }
    
    private double getNormalValue(double degrees, double center, double left,
                                  double right, double lower, double upper)
    {        
        double d;
        double radians = Math.toRadians(degrees);
        if (degrees >= 0.0 && degrees < 90.0) {
            d = (center-right)*Math.cos(radians)*Math.cos(radians) +
                (center-upper)*Math.sin(radians)*Math.sin(radians);
        } else if (degrees >= 90.0 && degrees < 180.0) {
            d = (center-left)*Math.cos(radians)*Math.cos(radians) +
                (center-upper)*Math.sin(radians)*Math.sin(radians);
        } else if (degrees >= 180.0 && degrees < 270.0) {
            d = (center-left)*Math.cos(radians)*Math.cos(radians) +
                (center-lower)*Math.sin(radians)*Math.sin(radians);
        } else if (degrees >= 270.0 && degrees < 360.0) {
            d = (center-right)*Math.cos(radians)*Math.cos(radians) +
                (center-lower)*Math.sin(radians)*Math.sin(radians);
        } else { return 0; }
        
        return d;
    }
    
    private double calculateStreamFunctionNormalDerivative(int i, int j)
    {
        Node node = field[i][j];
        double degrees = node.degreeOfSurfaceNormal;
        if (i - 1 < 0 || j - 1 < 0 || i + 1 > xNodes-1 || j + 1 > xNodes-1)
            return 0;
        double d = getNormalValue(degrees, node.psi,
                                  field[i-1][j].psi, field[i+1][j].psi,
                                  field[i][j-1].psi, field[i][j+1].psi);
        return 2*d/(cellSize*cellSize);
    }
    
    public void initializeStreamFunctionSolution()
    {
        double F = 1.4;
        applyGaussSeidelIteration(F);
    }
    
    public void advanceTimeStep()
    {
        double F = 1.4;
        advanceVorticity();
        applyGaussSeidelIteration(F);
    }
    
    public double[][] getVorticityData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                data[i][yNodes-j-1] = field[i][j].w;
            }
        }
        return data;
    }
    
    public double[][] getTemperatureData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                data[i][yNodes-j-1] = field[i][j].T;
            }
        }
        return data;
    }
    
    public double[][] getStreamFunctionData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                data[i][yNodes-j-1] = field[i][j].psi;
            }
        }
        return data;
    }
    
    public double[][] getXVelocityData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                data[i][yNodes-j-1] = field[i][j].u;
            }
        }
        return data;
    }
    
    public double[][] getYVelocityData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                data[i][yNodes-j-1] = field[i][j].v;
            }
        }
        return data;
    }
    
    public double[][] getNodeTypeData()
    {
        double[][] data = new double[xNodes][yNodes];
        for (int i = 0; i < xNodes; i++) {
            for (int j = 0; j < yNodes; j++) {
                switch (field[i][j].type) {
                    case STREAM:
                        data[i][yNodes-j-1] = 0.0;
                        break;
                    case OUTFLOW:
                    case INFLOW:
                    case FREELID:
                        data[i][yNodes-j-1] = 10.0;
                        break;
                    case SOLIDWALL:
                    case SOLIDSURFACE:
                        data[i][yNodes-j-1] = 20.0;
                        break;
                    case SOLIDINTERIOR:
                        data[i][yNodes-j-1] = 30.0;
                        break;
                    default:
                        break;
                }
            }
        }
        return data;
    }
    
    public static void main(String[] args)
    {
    }
}