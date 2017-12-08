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

shinyjs.hidePlotTooltip = function(){
     $("#ggvis-tooltip").hide();
};