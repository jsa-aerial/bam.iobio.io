// extending Thomas Down's original BAM js work

function Bam (bamUri, options, cb) {
    if (!(this instanceof arguments.callee)) {
        throw new Error("Constructor may not be called as a function");
    }

    this.bamUri = bamUri;
    this.options = options; // *** add options mapper ***
    // test if file or url
    if (typeof(this.bamUri) == "object") {
        var me = this;
        this.sourceType = "file";
        this.bamR = new readBinaryBAM(
          options.bai, bamUri,
          function(_) {
            me.bamR.bamFront(function(head, refs){
              me.bam = me.bamR;
              me.header = {sq: new2oldRefs(refs), head: head};
              if (cb) cb.call(me);
            });
          });
    } else if ( this.bamUri.slice(0,4) == "http" ) {
        this.sourceType = "url"
    }

    // set iobio servers
    this.iobio = {}
    this.iobio.bamtools = "ws://bamtools.iobio.io";
    this.iobio.samtools = "ws://samtools.iobio.io";
    this.iobio.bamReadDepther = "ws://bamReadDepther.iobio.io";
    this.iobio.bamMerger = "ws://bammerger.iobio.io";
    this.iobio.bamstatsAlive = "ws://bamstatsalive.iobio.io";
//      this.iobio.bamtools = "ws://localhost:8061";
//      this.iobio.samtools = "ws://localhost:8060";
//      this.iobio.bamReadDepther = "ws://localhost:8021";
//      this.iobio.bamMerger = "ws://localhost:8030";
//      this.iobio.bamstatsAlive = "ws://localhost:7100"
//      this.iobio.bamstatsAlive = "ws://localhost:7101"
    return this;
}


$.extend(Bam.prototype,{


   backstreamLclBamRegs: function (service, uid, refsNregions) {
       var allcnt = refsNregions.length;
       var regcnt = 0;
       var totcnt = 0;
       var bamR = this.bamR;
       backstream(
           service, uid,
           function(stream, opts, client){
               console.log(">>> Got STREAM event", opts)
               stream.write(bamR.headUba.buffer);
               bamR.throttledRegions2BAM(
                   refsNregions,
                   function(bgzfblks, contfn, regmap){
                       console.log("**** ", regcnt, bgzfblks)
                       if (bgzfblks && regcnt < allcnt) {
                           bgzfblks.forEach(
                               function(bgzfblk){stream.write(bgzfblk.buffer)});
                           totcnt = totcnt + bgzfblks.reduce(
                               function(S, x){return S+x.length},0);
                           regcnt = regcnt + 1;
                           contfn.call(this, regmap);
                       } else {
                           stream.write(EOFblk.buffer);
                           stream.end();
                           stream.destroy();
                           //client.close();
                           console.log("FINISHED")}})
           });
   },

   _getBamUrl: function(name, start, end) {
      return this._getBamRegionsUrl([ {'name':name,'start':start,'end':end} ]);
   },

   _getBamRegionsUrl: function(regions) {
      if ( this.sourceType == "url") {
         var regionStr = "";
         regions.forEach(function(region) { regionStr += " " + region.name + ":" + region.start + "-" + region.end });
         var url = this.iobio.samtools + "?cmd= view -b " + this.bamUri + regionStr + "&encoding=binary";

      } else {
         var me = this;
         var connectionID = makeuid();
         var url = "http://client?&id=" + connectionID;
         // setup backstream of binary sliced local bam
         me.backstreamLclBamRegs.call(
             me, me.iobio.bamstatsAlive,
             connectionID, regions);
      };
      return encodeURI(url);
   },

   _generateExomeBed: function(id) {
      var bed = "";
      var readDepth = this.readDepth[id];
      var start, end;
      var sum =0;
      // for (var i=0; i < readDepth.length; i++){
      //    sum += readDepth[i].depth;
      // }
      // console.log("avg = " + parseInt(sum / readDepth.length));
      // merge contiguous blocks into a single block and convert to bed format
      for( var i=0; i < readDepth.length; i++){
         if (readDepth[i].depth < 20) {
            if (start != undefined)
               bed += id + "\t" + start + "\t" + end + "\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
            start = undefined;
         }
         else {
            if (start == undefined) start = readDepth[i].pos;
            end = readDepth[i].pos + 16384;
         }
      }
      // add final record if data stopped on non-zero
      if (start != undefined)
         bed += id + "\t" + start + "\t" + end + "\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"
      return bed;
   },

   _tmp: function(ref, regions, bed){
      var me = this;
      var bedRegions = [];
      var a = this._bedToCoordinateArray(ref, bed);
      regions.forEach(function(reg) {
         var start = reg.start
         var length = reg.end - reg.start;
         var ci = me._getClosestValueIndex(a, reg.start); // update lo start value
         var maxci = a.length;
         while(length > 0 && ci < maxci) {
            var newStart,newEnd;

            // determine start position
            if ( a[ci].start <= start ) newStart = start;
            else newStart = a[ci].start;

            // determine end position
            if ( a[ci].end >=  newStart+length ) newEnd = newStart+length
            else { newEnd = a[ci].end; ci += 1; }

            // update length left to sample
            length -= newEnd - newStart;
            // push new regions
            bedRegions.push({ name:reg.name, 'start':newStart, 'end':newEnd});
         }
      })
      return bedRegions;
   },

   _mapToBedCoordinates: function(ref, regions, bed) {
      var a = this._bedToCoordinateArray(ref, bed);
      var a_i = 0;
      var bedRegions = [];
      if (a.length == 0) {
         alert("Bed file doesn't have coordinates for reference: " + regions[0].name + ". Sampling normally");
         return null;
      }
      regions.forEach(function(reg){
         for (a_i; a_i < a.length; a_i++) {
            if (a[a_i].end > reg.end)
               break;

            if (a[a_i].start >= reg.start)
               bedRegions.push( {name:reg.name, start:a[a_i].start, end:a[a_i].end})
         }
      })
      return bedRegions
   },

   _bedToCoordinateArray: function(ref, bed) {
      var me = this;
      var a = [];
      bed.split("\n").forEach(function(line){
        if (line[0] == '#' || line == "") return;

        var fields = line.split("\t");
        if (me._referenceMatchesBed(ref, fields[0])) {
           a.push({ chr:ref, start:parseInt(fields[1]), end:parseInt(fields[2]) });
        }
      });
      return a;
   },

   _referenceMatchesBed: function(ref, bedRef) {
      if (ref == bedRef) {
        return true;
      }
      // Try stripping chr from reference names and then comparing
      ref1 = ref.replace(/^chr?/,'');
      bedRef1 = bedRef.replace(/^chr?/,'');

      return (ref1 == bedRef1);
   },

   _getClosestValueIndex: function(a, x) {
       var lo = -1, hi = a.length;
       while (hi - lo > 1) {
           var mid = Math.round((lo + hi)/2);
           if (a[mid].start <= x) {
               lo = mid;
           } else {
               hi = mid;
           }
       }
       if (lo == -1 ) return 0;
       if ( a[lo].end > x )
           return lo;
       else
           return hi;
   },

   getReferencesWithReads: function(callback) {
      var me = this;
      if (this.sourceType == 'url') {

      } else {
          refsWithReads.call(me, callback, true);
      }
   },


   // *** bamtools functionality ***

   convert: function(format, name, start, end, callback, options) {
      // Converts between BAM and a number of other formats
      if (!format || !name || !start || !end)
         return "Error: must supply format, sequenceid, start nucleotide and end nucleotide"

      if (format.toLowerCase() != "sam")
         return "Error: format + " + options.format + " is not supported"
      var me = this;
      this.fetch(name, start, end, function(data,e) {
         if(options && options.noHeader)
            callback(data, e);
         else {
            me.getHeader(function(h) {
               callback(h.toStr + data, e);
            })
         }
      }, { 'format': format })
   },

   count: function() {
      // Prints number of alignments in BAM file(s)
   },

   coverage: function() {
      // Prints coverage statistics from the input BAM file
   },

   filter: function() {
      // Filters BAM file(s) by user-specified criteria
   },


   estimateBaiReadDepth: function(callback) {
      var me = this, readDepth = {};
      me.readDepth = {};

      function cb() {
         if (me.header) {
            for (var id in readDepth) {
              if (readDepth.hasOwnProperty(id))
              var name = me.header.sq[parseInt(id)].name;
               if ( me.readDepth[ name ] == undefined){
                  me.readDepth[ name ] = readDepth[id];
                  callback( name, readDepth[id] );
               }
            }
         }
      }

      me.getHeader(function(header) {
         if (Object.keys(me.readDepth).length > 0)
            cb();
      });

      if ( Object.keys(me.readDepth).length > 0 )
         callback(me.readDepth)
      else if (me.sourceType == 'url') {
         var currentSequence;
         doService(
             me.iobio.bamReadDepther,
             encodeURI( me.iobio.bamReadDepther + '?cmd=-i ' + me.bamUri + ".bai"),
             function(data) {
                if (data == undefined) {
                    cb();
                } else {
                    data = data.split("\n");
                    for (var i=0; i < data.length; i++)  {
                        if ( data[i][0] == '#' ) {
                            if ( Object.keys(readDepth).length > 0 ) { cb() };
                            currentSequence = data[i].substr(1);
                            readDepth[currentSequence] = [];
                        } else if (data[i] != "") {
                            var d = data[i].split("\t");
                            readDepth[currentSequence].push({ pos:parseInt(d[0]), depth:parseInt(d[1]) });
                        };
                    };
                };
             });

      } else if (me.sourceType == "file") {
          var bamR = me.bamR;
          var refs = bamR.refsWithReads().map(function(x){return x[1].name;});
          refs.forEach(function(r){
              me.readDepth[r] = mapSegCoverage(
                  bamR.baiR, bamR.refName2Index(r),
                  function(x){
                    return (x.depth > 0) ? x : undefined;
                  }).sort(function(l,r){return l.pos - r.pos});
              callback(r, me.readDepth[r]);
          });
      }

   },


   getHeader: function(callback) {
      var me = this;
      if (me.header)
         callback(me.header);
      else if (me.sourceType == "file")
          me.bamR.bamFront(function(head, refs) {
              me.header = {sq: new2oldRefs(refs), head: head};
              callback.call(me, me.header);
          });
      else {
         var url = encodeURI( me.iobio.samtools + '?cmd=view -H ' + this.bamUri)
         var rawHeader = ""
         doService(
             me.iobio.samtools, url,
             function(data) {
               if (data){
                 rawHeader += data;
               } else {
                   me.setHeader(rawHeader);
                   callback( me.header);
               };
             });
      };
   },

   setHeader: function(headerStr) {
      var header = { sq:[], toStr : headerStr };
      var lines = headerStr.split("\n");
      for ( var i=0; i<lines.length > 0; i++) {
         var fields = lines[i].split("\t");
         if (fields[0] == "@SQ") {
            var fHash = {};
            fields.forEach(function(field) {
              var values = field.split(':');
              fHash[ values[0] ] = values[1]
            })
            header.sq.push({name:fHash["SN"], end:1+parseInt(fHash["LN"])});
         }
      }
      this.header = header;
   },

   index: function() {
      // Generates index for BAM file
   },

   merge: function() {
      // Merge multiple BAM files into single file
   },

   random: function() {
      // Select random alignments from existing BAM file(s), intended more as a testing tool.
   },

   resolve: function() {
      // Resolves paired-end reads (marking the IsProperPair flag as needed)
   },

   revert: function() {
      // Removes duplicate marks and restores original base qualities
   },

   sort: function() {
      // Sorts the BAM file according to some criteria
   },

   split: function() {
      // Splits a BAM file on user-specified property, creating a new BAM output file for each value found
   },



   getBedRegions: function (ref, regions, options) {
       var bed = me._generateExomeBed(options.sequenceNames[0]);
       var bedRegions = undefined;
       // map random region coordinates to bed coordinates
       if (bed != undefined) {
           bedRegions =
               me._mapToBedCoordinates(ref.name, regions, bed);
       };
       return {bed: bed, bedRegions: bedRegions};
   },

   getSamplingRegions: function(refs, options) {
       var me = this;
       var regions = samplingRegions(refs, options);
       var regionInfo = {regions: regions};
       if (options.exomeSampling) {
           var bedInfo = me.getBedRegions(refs[0], regions, options);
           regionInfo.bedRegions = bedInfo.bedRegions;
           regionInfo.bed = bedInfo.bed;
           options.bed = bedInfo.bed; // Backward compatible...
       }
       return regionInfo;
   },

    _getBamUrl: function(name, start, end) {
        return this._getBamRegionsUrl([ {'name':name,'start':start,'end':end} ]);
    },



   stats: function(name, start, end, callback) {
      // Prints some basic statistics from input BAM file(s)
      var url =
           encodeURI(this.iobio.bamstatsAlive + '?cmd=-u 1000 -s ' + start +
                     " -l " + parseInt(end-start) + " " +
                     encodeURIComponent(this._getBamUrl(name,start,end)));
      var buffer = "";
      doService(
          this.iobio.bamstatsAlive, url,
          function(data) {
            if (data == undefined) return;
            var success = true;
            try {
              var obj = JSON.parse(buffer + data)
            } catch(e) {
              success = false;
              buffer += data;
            }
            if(success) {
              buffer = "";
              callback(obj);
            }
          });
   },


   bamStatsAliveSamplingURL: function (SQs, options) {
       var me = this;
       var regionInfo = me.getSamplingRegions(SQs, options);
       var regions = regionInfo.regions;
       var bedRegions = regionInfo.bedRegions;
       var regStr = JSON.stringify((bedRegions || regions).map(
           function(r) {
               return {start:r.start, end:r.end, chr:r.name};
           }));
       var urlOpts =
           (me.sourceType == "url") ?
               '?cmd=-u 500 -r \'' + regStr + '\' ' :
               "?protocol=websocket&cmd=-u 500 ";

       return encodeURI(me.iobio.bamstatsAlive + urlOpts +
                        encodeURIComponent(me._getBamRegionsUrl(regions)));
   },

   sampleStats: function(callback, options) {
      // Prints some basic statistics from sampled input BAM file(s)
      options = $.extend({
         binSize : 40000, // defaults
         binNumber : 20,
         start : 1,
      },options);
      var me = this;

      function goSampling(SQs) {
         var url = me.bamStatsAliveSamplingURL(SQs, options);

         var buffer = "";
         doService(
             me.iobio.bamstatsAlive, url,
             function(datas) {
               datas.split(';').forEach(function(data) {
                 if (data == undefined) {
                   if (options.onEnd != undefined) options.onEnd();
                   return;
                 } else {
                   var success = true;
                   try {
                     var obj = JSON.parse(buffer + data)
                   } catch(e) {
                     success = false;
                     buffer += data;
                   };
                   if(success) {
                     buffer = "";
                     callback(obj);
                   };
                 };
               });
             });
      }

      if ( options.sequenceNames != undefined && options.sequenceNames.length == 1 && options.end != undefined && me.sourceType != "file") {
         goSampling([{name:options.sequenceNames[0], end:options.end}]);
      } else if (options.sequenceNames != undefined && options.sequenceNames.length == 1){
         this.getHeader(function(header){
            var sq;
            $(header.sq).each( function(i,d) {
               if(d.name == options.sequenceNames[0])
               sq = d;
            })
            goSampling([sq]);
         });
      } else if (me.sourceType == "file") {
          goSampling(me.bamR.refsWithReads().map(function(x) {return x[1]}));
      } else {
         this.getHeader(function(header){
            goSampling(header.sq);
         });
      }
   }

});
