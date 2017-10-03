var rtermFormat = require('./rtermFormat');



const events = {
    update: function(ui) {
//        randomTerms(ui, this);
    },

    onChange_factors: function(ui) {
        randomTerms(ui, this);
        plotTerms(ui, this);
    },

    onChange_covariates: function(ui) {
        randomTerms(ui, this);
//        plotTerms(ui, this);

    },
    onChange_cluster: function(ui) {
        randomTerms(ui, this);
    },

    onChange_possibleTerms: function(ui) {
     console.log("possible changed");
    },

    onChange_randomTerms: function(ui) {
        filterRandomTerms(ui, this);
    },
    onChange_fixedTerms: function(ui) {
      filterFixedTerms(ui, this);
    },


    onEvent_randomTerms_preprocess: function(ui, data) {
        for(var j = 0; j < data.items.length; j++) {
          data.items[j].value.raw=data.items[j].value.toString();
        }
}
};

var randomTerms = function(ui, context) {

    var factorList = context.cloneArray(ui.factors.value(), []);
    var covariatesList = context.cloneArray(ui.cov.value(), []);

    var variableList = factorList.concat(covariatesList);
    var termsList=[];
    termsList = context.getCombinations(variableList);
    context.sortArraysByLength(termsList);
    ui.fixedSupplier.setValue(context.valuesToItems(termsList, FormatDef.term));
    var clusterList = context.cloneArray(ui.cluster.value(), []);
    if (clusterList.length<1) {
                ui.modelSupplier.setValue(context.valuesToItems([], rtermFormat));   
                return;
    }
    termsList.unshift(["Intercept"]);
    var alist=[];
    for (var i=0; i < clusterList.length; i++) {
     for (var j = 0; j < termsList.length; j++) {
       var item=context.cloneArray(termsList[j]);
       item[item.length]=clusterList[i];
       alist.push(item);
//       console.log(item);
     }
    }
    context.sortArraysByLength(alist);
    var formatted=context.valuesToItems(alist, rtermFormat);

    ui.modelSupplier.setValue(formatted);
};


var filterRandomTerms = function(ui, context) {
    var termsList = context.cloneArray(ui.randomTerms.value(), []);
    var unique = termsList.filter((v, i, a) => a.indexOf(v) === i); 
    if (unique.length!=termsList.length)
      ui.randomTerms.setValue(unique);

};


var filterFixedTerms = function(ui, context) {
    var termsList = context.cloneArray(ui.fixedTerms.value(), []);
    var oList=context.valuesToItems(termsList, FormatDef.term);
    var aList= listToStrings(oList);
    var unique = aList.filter((v, i, a) => a.indexOf(v) === i); 
    var delta=aList.length-unique.length;
    var pos=aList.map( function(item) {return aList.indexOf(item);});
        pos = pos.filter((v, i, a) => a.indexOf(v) === i);
    var newList=[];    
        pos.forEach(function(val) {newList.push(termsList[val])});
        
    if (newList.length!=termsList.length) {
          ui.fixedTerms.setValue(newList);
    }
    };
    

var listToStrings = function(obj) {
    var aList = [];
       for (var i=0; i < obj.length; i++) {
          aList.push(obj[i].value.toString());
       }
       return(aList);
       
    };




module.exports = events;
