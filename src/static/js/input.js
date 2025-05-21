window.routeInput = 
{
    start: { name: "鸡鸣寺站5号口", type: "STA", feature_type: "Point"},
    end: { name: "解放门", type: "TOU", feature_type: "Point"},
    segments: {
      segment_1: {
        segment_start: { name: "鸡鸣寺站5号口", type: "STA", feature_type: "Point"},
        segment_end: { name: "鸡鸣寺", type: "TOU", feature_type: "Polygon"},
        path_instructions: [
          { action:"out"},
          { action:"Turn",direction: "right"},
          { transport:"walk",distance: "500m" }
        ]
      },
      segment_2: {
          segment_start: { name: "鸡鸣寺", type: "TOU", feature_type: "Polygon"},
          segment_end: { name: "解放门", type: "TOU", feature_type: "Point"},
          path_instructions: [
            { transport:"walk",distance: "100m" }
          ]
        }
    }
  };
