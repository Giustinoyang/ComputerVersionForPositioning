% Demo for Edge Boxes (please see readme.txt first).

%% load pre-trained edge detection model and set opts (see edgesDemo.m)
model=load('models/forest/modelBsds'); model=model.model;
model.opts.multiscale=0; model.opts.sharpen=2; model.opts.nThreads=4;

%% set up opts for edgeBoxes (see edgeBoxes.m)
opts = edgeBoxes;
opts.alpha = .65;     % step size of sliding window search
opts.beta  = .75;     % nms threshold for object proposals
opts.minScore = .01;  % min score of boxes to detect
opts.maxBoxes = 1e4;  % max number of boxes to detect

%% detect Edge Box bounding box proposals (see edgeBoxes.m)
test_path = '~/Documents/sandbox/test_edge_boxes/test/';
d = dir([test_path, '/*.jpg']);
for i = 1 : numel(d)
    I = imread([test_path, '/', d(i).name]);
    bbs = edgeBoxes(I,model,opts);

    n = 25;
    bbtmp = [bbs(1:n, :), ones(n, 1)];
    clf; hs = bbGt('showRes',I,[],bbtmp);

    colors = 'ymcrgbwk';

    for j = 1:n
        set(hs(2*j-1), 'EdgeColor', colors(mod(j,numel(colors)) + 1));
    end
    
    F = getframe;
    imwrite(F.cdata, [test_path, '/../results_matlab/', d(i).name]);
end