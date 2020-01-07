import io.jbotsim.ui.JViewer;

public class GeneralColo {
    public static void main(String[] args){
        MyTopology tp = new MyTopology();
        tp.setDefaultNodeModel(GeneralNode.class);
        new JViewer(tp);
    }
}