/**********************************************
  * Ryan Shea
  * MAE 423 - Heat Transfer
  * Final Project
  **********************************************/

import java.awt.Color;

public class FlowVisualizer
{
    private static void drawSite(double value, int i, int j, double cellSize,
                                 double minValue, double maxValue,
                                 double height, double offset, double colorRangePercentage)
    {
        if (value > maxValue) value = maxValue;
        if (value < minValue) value = minValue;
        Color color = ColorHelper.numberToColor((value-minValue)/(maxValue-minValue)*colorRangePercentage);
        StdDraw.setPenColor(color);
        double x = cellSize*i - cellSize/2;
        double y = height - (cellSize*j - cellSize/2) + offset;
        StdDraw.filledSquare(x, y, cellSize/2);
    }
    
    public static void main(String[] args)
    {
        double scale = 1000.0;
        double width = 500.0/scale;
        double height = 200.0/scale;
        double cellSize = 1.0/scale;
        double screenWidth = width/cellSize;
        double screenHeight = height/cellSize;
        double screenScaleFactor = 1.5*scale;
        double cylinderDiameter = height/4.0;
        
        StdDraw.setCanvasSize((int) (width*screenScaleFactor),
                              (int) (height*screenScaleFactor*2));
        StdDraw.setXscale(0.0, width);
        StdDraw.setYscale(0.0, height*2);
        StdDraw.show(0);
        
        double flowTemperature = 300.0;
        double surfaceTemperature = 400.0;
        double reynoldsNumber = 200.0;
        double kinematicViscosity = 0.000001569;
        double alpha = 0.0002216;
        Flow field = new Flow(width, height, cellSize,
                              flowTemperature, reynoldsNumber,
                              kinematicViscosity, alpha, cylinderDiameter);
        field.addSolidCylinder(width*1/4, height/2.0, cylinderDiameter,
                               surfaceTemperature);
        /*field.addSolidBox(width*1/6, height*0.25, width/40.0, height/10.0,
                          surfaceTemperature);
        field.addSolidBox(width*1/6, height*0.75, width/40.0, height/10.0,
                          surfaceTemperature);
        field.addSolidCylinder(width*2/6, height*1/2, height/10.0,
                               surfaceTemperature);
        field.addSolidBox(width*3/6, height*0.25, width/40.0, height/10.0,
                          surfaceTemperature);
        field.addSolidBox(width*3/6, height*0.75, width/40.0, height/10.0,
                          surfaceTemperature);
        field.addSolidCylinder(width*4/6, height*1/2, height/10.0,
                               surfaceTemperature);
        field.addSolidBox(width*5/6, height*0.25, width/40.0, height/10.0,
                          surfaceTemperature);
        field.addSolidBox(width*5/6, height*0.75, width/40.0, height/10.0,
                          surfaceTemperature);*/
        field.initializeStreamFunctionSolution();
        
        //double[][] data1 = field.getTemperatureData();
        //double[][] data1 = field.getXVelocityData();
        double[][] data1 = field.getStreamFunctionData();
        int xNodes1 = data1.length;
        int yNodes1 = data1[0].length;
        double maxValue1 = 0;
        for (int i = 0; i < xNodes1; i++) {
            for (int j = 0; j < yNodes1; j++) {
                if (maxValue1 < data1[i][j])
                    maxValue1 = data1[i][j];
            }
        }
        
        //double[][] data2 = field.getStreamFunctionData();
        double[][] data2 = field.getVorticityData();
        //double[][] data2 = field.getXVelocityData();
        int xNodes2 = data2.length;
        int yNodes2 = data2[0].length;
        double maxValue2 = Double.NEGATIVE_INFINITY;
        double minValue2 = Double.POSITIVE_INFINITY;
        
        int iterations = 0;
        double time = 0;
        while (true)
        {            
            for (int i = 0; i < xNodes2; i++) {
                for (int j = 0; j < yNodes2; j++) {
                    if (maxValue2 < data2[i][j])
                        maxValue2 = data2[i][j];
                    if (minValue2 > data2[i][j])
                        minValue2 = data2[i][j];
                }
            }
            
            for (int i = 0; i < xNodes1; i++) {
                for (int j = 0; j < yNodes1; j++) {
                    FlowVisualizer.drawSite(data1[i][j], i, j, cellSize,
                                            0, maxValue1, height, 0, 100);
                }
            }
            
            for (int i = 0; i < xNodes2; i++) {
                for (int j = 0; j < yNodes2; j++) {
                    FlowVisualizer.drawSite(data2[i][j], i, j, cellSize,
                                            minValue2/10, maxValue2/10, height, height, 50);
                }
            }
            StdDraw.setPenColor(Color.BLACK);
            StdDraw.filledRectangle(width/2, height, height/5, height/20);
            StdDraw.setPenColor(Color.WHITE);
            time = iterations*field.dt;
            //StdDraw.text(width/2, height, String.format("%1$,.2f", time));
            StdDraw.text(width/2, height, Integer.toString(iterations));
            StdDraw.show(50);
            if (iterations >= 15000) break;
            for (int i = 0; i < 50; i++) {
                field.advanceTimeStep();
                iterations++;
            }
            //data1 = field.getStreamFunctionData();
            data2 = field.getVorticityData();
            //data1 = field.getXVelocityData();
            data1 = field.getStreamFunctionData();
            //data1 = field.getTemperatureData();
        }
    }
}