shinyjs.coinPlot = function(){
    $.each($('#g64').children(),function(){$(this).attr('d','M0,4.675581A6.6755811781245455,1.6755811781245455 0 1,1 0,-3.675581A6.6755811781245455,1.6755811781245455 0 1,1 0,4.6755811781245455Z')});
};


shinyjs.changeTree = function(params){
    eval("$('#tree').jstree(true).settings.core.data="+params);
    $('#tree').jstree(true).refresh();
};

shinyjs.open = function(){
    $('#tree').jstree(true).open_all();
};

shinyjs.deselect = function(){
    $('#tree').jstree(true).deselect_all();
};

shinyjs.setDefaultTree = function(){
    $.jstree.defaults.checkbox.three_state= false;
};

shinyjs.openStaticTree = function(){
    $('#staticRegionTree').jstree(true).open_all();
};