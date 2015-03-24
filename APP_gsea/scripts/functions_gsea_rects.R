rects.center = function(rectangle) {
    r = rectangle
    x = (r$xmax + r$xmin)/2
    y = (r$ymax + r$ymin)/2
    return(c(x, y))
}

rects.draw = function(rectangles, ...) {
    # draw a single rectangle from rects.get or list of rectangles from rects.get extra args will be passed to
    # rect()
    if(!is.list(rectangles)){
      stop('input must be a list! either a list with names c("xmin", "xmax", "ymin", "ymax") or a list of such lists.')
    }
    rects = rectangles
    if (!is.null(names(rects)) & length(names(rects)) == 4) {
      if(all(names(rects) == c("xmin", "xmax", "ymin", "ymax")))
        rects = list(rectangles)
    }
    for (i in 1:length(rects)) {
        r = rects[[i]]
        rect(xleft = r$xmin, ybottom = r$ymin, xright = r$xmax, ytop = r$ymax, ...)
    }
}

rects.get = function(vertLines = log2(c(1/2, 2)), horLines = 1, XMIN = -10, XMAX = 10, 
    YMIN = 0, YMAX = 6) {
    # with no args, returns a predefined set of rectangles regions in x are from XMIN to XMAX split by vertLines
    # regions in y are from YMIN to YMAX split by horLines XMIN = -4 XMAX = 4 YMIN = 0 YMAX = 6 vertLines =
    # c(1/4, 1/2.5, 1/1.5, 1.5, 2.5,4) vertLines = log2(vertLines) horLines = 2 horLines = log2(horLines)
    
    all_rects = list()
    add_rect = function(all_rects, xmin, xmax, ymin, ymax) {
        all_rects[[length(all_rects) + 1]] = list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
        return(all_rects)
    }
    
    add_rowRects = function(all_rects, ymin, ymax) {
        all_rects = add_rect(all_rects, xmin = XMIN, xmax = vertLines[1], ymin = ymin, ymax = ymax)
        for (v in 1:(length(vertLines) - 1)) {
            all_rects = add_rect(all_rects, xmin = vertLines[v], xmax = vertLines[v + 1], ymin = ymin, ymax = ymax)
        }
        all_rects = add_rect(all_rects, xmin = vertLines[length(vertLines)], xmax = XMAX, ymin = ymin, ymax = ymax)
        return(all_rects)
    }
    all_rects = add_rowRects(all_rects, horLines[1], YMAX)
    all_rects = add_rowRects(all_rects, YMIN, horLines[1])
    return(all_rects)
}

rects.countIn = function(x, y, rectangle) {
    if (length(x) != length(y)) {
        print("data length mismatch!")
        stop
    }
    xmin = rectangle$xmin
    xmax = rectangle$xmax
    ymin = rectangle$ymin
    ymax = rectangle$ymax
    # x = get_x(dat, end, start) y = get_y(dat, end, start)
    is_inRect = y <= ymax & y >= ymin & x <= xmax & x >= xmin
    # print(sum(is_inRect)/length(x))
    return(is_inRect)
} 

