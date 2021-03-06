import io.jbotsim.core.Topology;
import io.jbotsim.core.Node;

import java.util.List;

public class MyTopology extends Topology {

    @Override
    public void start() {
        List<Node> nodeList = getNodes();
        MyNode.nb_nodes = nodeList.size();
        for (Node node: nodeList ) {
            GeneralNode myNode = (GeneralNode) node;
            myNode.init_delta();
        }
        super.start();
    }
}
